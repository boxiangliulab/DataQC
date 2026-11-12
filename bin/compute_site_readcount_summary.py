#!/usr/bin/env python3
"""
compute_site_readcount_summary.py

Splicing-specific DataQC summary that the expression QC has no analogue for: per-splice-site
RAW read-count statistics (missingness = zero-count fraction, plus total/mean/median/min/max
read DEPTH on the denominator, and `std_ratio` = SD of the PSI RATIO), computed from the Stage-0
ISSAC `.site` files and restricted to the sites that survived QC (i.e. that are present in the
QC'd combined splice matrix).

`std_ratio` is the standard deviation of the per-sample splice ratio `(numerator + 0.5) /
(denominator + 0.5)` over the OBSERVED samples (denominator >= 1), matching the smoothed PSI that
Stage-0 writes into the phenotype matrix. Count columns remain based on the denominator (depth).

Mirrors issac_bulk.nf COMPUTE_SUMMARY_STATS / COMPUTE_SUMMARY_STATS_IR, but runs inside 1_DataQC so
the stats are tied to the analysed site set. Streams the (large) `.site` files line by line
(memory-light). Junction and IR are summarised separately.

`.site` format: line 1 = space-separated sample-ID header (no ':'), each subsequent line =
  <site_id> <include:total> <include:total> ...   (site_id like chr:strand:pos)
The QC'd feature IDs (keep_sites) look like `junction_<site_id>` / `IR_<site_id>`; we strip the
`junction_`/`IR_` prefix to join to the `.site` site_id.
"""

import argparse
import os
import sys
import numpy as np


def load_keep_sites(keep_sites_path):
    """Read QC'd feature IDs (junction_<id> / IR_<id>); return two sets of bare site_ids."""
    junc, ir = set(), set()
    if not keep_sites_path or not os.path.exists(keep_sites_path):
        return None, None  # None => keep all sites
    with open(keep_sites_path) as f:
        for line in f:
            fid = line.strip().split('\t')[0].strip()
            if not fid or fid.lower() in ('feature', 'id', 'junction', 'gene'):
                continue
            if fid.startswith('junction_'):
                junc.add(fid[len('junction_'):])
            elif fid.startswith('IR_'):
                ir.add(fid[len('IR_'):])
    return (junc or None), (ir or None)


def summarise_site_file(site_path, keep_set, modality, out_tsv):
    """Stream a .site file; write per-site stats TSV; return the list of missingness/depth arrays."""
    miss_all, depth_all = [], []
    n_sites = 0
    with open(site_path) as fin, open(out_tsv, 'w') as fout:
        fout.write('\t'.join(['site_id', 'missingness', 'total_counts', 'mean_count',
                              'std_ratio', 'median_count', 'min_count', 'max_count',
                              'modality']) + '\n')
        for line in fin:
            fields = line.rstrip('\n').split(' ')
            if len(fields) < 2:
                continue
            site_id = fields[0]
            totals = []       # denominator (read depth) per sample
            ratios = []       # smoothed PSI ratio over OBSERVED samples (denominator >= 1)
            for fld in fields[1:]:
                if ':' in fld:
                    parts = fld.split(':')
                    try:
                        num = int(parts[0]); tot = int(parts[1])
                    except (ValueError, IndexError):
                        num = 0; tot = 0
                    totals.append(tot)
                    if tot >= 1:
                        ratios.append((num + 0.5) / (tot + 0.5))
            if not totals:
                continue  # header line (sample IDs) / malformed -> skip
            if keep_set is not None and site_id not in keep_set:
                continue
            c = np.array(totals, dtype=float)
            r = np.array(ratios, dtype=float)
            n = len(c)
            miss = float(np.sum(c == 0)) / n
            std_ratio = float(r.std(ddof=1)) if r.size > 1 else 0.0
            fout.write('\t'.join([
                site_id, f"{miss:.6f}", str(int(c.sum())), f"{c.mean():.6f}",
                f"{std_ratio:.6f}", f"{np.median(c):.6f}",
                f"{c.min():.6f}", f"{c.max():.6f}", modality]) + '\n')
            miss_all.append(miss)
            depth_all.append(c.mean())
            n_sites += 1
    return n_sites, np.array(miss_all), np.array(depth_all)


def agg_line(modality, n_sites, miss, depth):
    def q(a, p):
        return float(np.percentile(a, p)) if a.size else float('nan')
    return [modality, str(n_sites),
            f"{q(miss,50):.4f}", f"{q(miss,25):.4f}", f"{q(miss,75):.4f}",
            f"{q(depth,50):.4f}", f"{q(depth,25):.4f}", f"{q(depth,75):.4f}"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--junction_site', default=None, help='Stage-0 03_phenotype/<cohort>.site')
    ap.add_argument('--ir_site', default=None, help='Stage-0 04_IR/<cohort>_IR.site')
    ap.add_argument('--keep_sites', default=None,
                    help='QC\'d feature list (junction_list.txt); if omitted, all sites are summarised')
    ap.add_argument('--output_dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    junc_keep, ir_keep = load_keep_sites(args.keep_sites)

    agg = [['modality', 'n_sites', 'missingness_median', 'missingness_q25', 'missingness_q75',
            'mean_depth_median', 'mean_depth_q25', 'mean_depth_q75']]

    if args.junction_site and os.path.exists(args.junction_site):
        n, miss, depth = summarise_site_file(
            args.junction_site, junc_keep, 'junction',
            os.path.join(args.output_dir, 'site_readcount_summary_junction.tsv'))
        agg.append(agg_line('junction', n, miss, depth))
        print(f"junction: summarised {n} sites", file=sys.stderr)
    else:
        print("WARNING: junction .site not found; skipping junction read-count summary", file=sys.stderr)

    if args.ir_site and os.path.exists(args.ir_site):
        n, miss, depth = summarise_site_file(
            args.ir_site, ir_keep, 'IR',
            os.path.join(args.output_dir, 'site_readcount_summary_IR.tsv'))
        agg.append(agg_line('IR', n, miss, depth))
        print(f"IR: summarised {n} sites", file=sys.stderr)
    else:
        print("WARNING: IR .site not found; skipping IR read-count summary", file=sys.stderr)

    with open(os.path.join(args.output_dir, 'site_readcount_summary_overview.tsv'), 'w') as f:
        for row in agg:
            f.write('\t'.join(row) + '\n')
    print("Site read-count summary written to " + args.output_dir, file=sys.stderr)


if __name__ == '__main__':
    main()
