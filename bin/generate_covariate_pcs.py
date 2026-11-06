#!/usr/bin/env python3
"""
Generate CovariatePCs.txt for sQTL analysis
============================================
Merges genotype PCs, splicing PCs, and optionally GenSex into eQTLGen-compatible format.

Input files:
- GenotypePCs.txt: Tab-delimited, samples in rows, 10 PCs in columns (from DataQC)
- Splicing PCs file: Tab-delimited, samples in rows, PCs in columns (from ProcessSplicing.R)
- Sexcheck file (optional): Tab-delimited, contains STATUS column for GenSex

Output:
- CovariatePCs.txt: SampleID | GenPC1-10 | SplicingPC1-50 | GenSex
"""

import argparse
import os
import pandas as pd
import sys

def read_genotype_pcs(filepath, num_pcs=10):
    """
    Read genotype PCs from DataQC output.
    Format: Tab-delimited, first column is sample ID, followed by PC1-PC10
    """
    print(f"Reading genotype PCs from: {filepath}")
    df = pd.read_csv(filepath, sep='\t', dtype={0: str})

    # Get first column name (should be empty or sample ID column)
    first_col = df.columns[0]

    # Select first column + num_pcs columns
    df = df.iloc[:, :num_pcs+1]

    # Rename columns
    df.columns = ['SampleID'] + [f'GenPC{i}' for i in range(1, num_pcs+1)]

    print(f"  Loaded {len(df)} samples with {num_pcs} genotype PCs")
    return df

def read_splicing_pcs(filepath, num_pcs=50, prefix='SplicingPC'):
    """
    Read phenotype PCs (splicing or expression) from Process{Splicing,Expression}.R output.
    Format: Tab-delimited, samples in rows, PCs in columns; first column = sample ID (genotype-ID space),
    followed by PC1, PC2, ...  Columns are renamed to <prefix>1..<prefix>N (e.g. SplicingPC / ExpressionPC).
    """
    print(f"Reading phenotype PCs ({prefix}) from: {filepath}")
    df = pd.read_csv(filepath, sep='\t', dtype={0: str})

    # Check format - should be samples x PCs
    print(f"  Original format: {df.shape[0]} samples x {df.shape[1]-1} PCs")

    # Get first column name (should be 'sample' / 'Sample' / 'SampleID')
    first_col = df.columns[0]

    # Select first column + num_pcs columns
    num_available_pcs = min(num_pcs, df.shape[1] - 1)
    df = df.iloc[:, :num_available_pcs+1]

    # Rename columns to standard format
    df.columns = ['SampleID'] + [f'{prefix}{i}' for i in range(1, num_available_pcs+1)]

    print(f"  Loaded {len(df)} samples with {num_available_pcs} {prefix} PCs")
    return df

def read_gen_sex(filepath):
    """
    Read GenSex from SexCheck.txt file (from DataQC).
    Format: Tab-delimited, FID IID PEDSEX SNPSEX STATUS F ...
    Extract IID and SNPSEX columns (SNPSEX: 0=UNKNOWN, 1=MALE, 2=FEMALE)
    """
    print(f"Reading GenSex from: {filepath}")
    df = pd.read_csv(filepath, sep='\t', dtype={0: str, 1: str})

    # Check if SNPSEX column exists and is not all NA
    if 'SNPSEX' not in df.columns or df['SNPSEX'].isna().all():
        print("  WARNING: No valid SNPSEX column found, skipping GenSex")
        return None

    # Get IID (second column) and SNPSEX (4th column)
    # Column names should be: FID IID PEDSEX SNPSEX STATUS F ...
    if len(df.columns) < 4:
        print("  WARNING: SexCheck file has fewer than 4 columns, skipping GenSex")
        return None

    # Extract IID (column 2) and SNPSEX (column 4)
    # Using column indices: 1 for IID, 3 for SNPSEX (0-indexed)
    sex_df = df.iloc[:, [1, 3]].copy()
    sex_df.columns = ['SampleID', 'GenSex']

    # Convert GenSex to integer
    sex_df['GenSex'] = pd.to_numeric(sex_df['GenSex'], errors='coerce')

    print(f"  Loaded GenSex for {len(sex_df)} samples")
    print(f"    UNKNOWN (0): {(sex_df['GenSex'] == 0).sum()}")
    print(f"    MALE (1): {(sex_df['GenSex'] == 1).sum()}")
    print(f"    FEMALE (2): {(sex_df['GenSex'] == 2).sum()}")
    return sex_df

def read_additional_covariates(filepath):
    """
    Port of the reference MakeCovariateFile.R additional-covariate handling.
    - First column must identify the sample (renamed to 'SampleID').
    - Rows with any missing value are dropped (complete.cases).
    - Binary text covariates (exactly 2 levels) -> 0/1 (sorted-level factor, matching
      R as.numeric(as.factor(x)) - 1).
    - Multi-level text covariates -> one-hot dummy columns named <col>_<level> (all levels kept).
    - Numeric covariates kept as-is.
    """
    print(f"Reading additional covariates from: {filepath}")
    df = pd.read_csv(filepath, sep='\t', dtype={0: str})
    if df.columns[0] != 'SampleID':
        df = df.rename(columns={df.columns[0]: 'SampleID'})
    df['SampleID'] = df['SampleID'].astype(str)

    # complete.cases: drop any row with a missing value
    df = df.dropna(axis=0, how='any').reset_index(drop=True)

    out = df[['SampleID']].copy()
    for col in df.columns:
        if col == 'SampleID':
            continue
        s = df[col]
        if pd.api.types.is_numeric_dtype(s):
            out[col] = s.values
        elif s.nunique() == 2:
            # sorted-level factor -> 0/1
            out[col] = pd.Categorical(s).codes
        else:
            dummies = pd.get_dummies(s, prefix=col, prefix_sep='_').astype(int)
            out = pd.concat([out, dummies], axis=1)

    print(f"  Loaded {len(out)} samples, {out.shape[1] - 1} covariate column(s) after encoding")
    return out


def main():
    parser = argparse.ArgumentParser(
        description='Generate CovariatePCs.txt for sQTL analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--genotype_pcs', required=True,
                        help='Path to GenotypePCs.txt from DataQC')
    parser.add_argument('--splicing_pcs', required=True,
                        help='Path to splicing PCs file from Step 1')
    parser.add_argument('--gen_sex', default=None,
                        help='Path to sexcheck file from DataQC (optional)')
    parser.add_argument('--output', required=True,
                        help='Output CovariatePCs.txt file')
    parser.add_argument('--num_gen_pcs', type=int, default=10,
                        help='Number of genotype PCs to use (default: 10)')
    parser.add_argument('--pc_prefix', default='SplicingPC',
                        help="Column-name prefix for the phenotype PCs (e.g. SplicingPC or ExpressionPC)")
    parser.add_argument('--num_splicing_pcs', type=int, default=50,
                        help='Number of splicing PCs to use (default: 50)')
    parser.add_argument('--additional_covariates', default=None,
                        help='Optional cohort-specific covariate file (tab-delimited, first column '
                             'SampleID). Categorical text covariates are auto-encoded (binary -> 0/1, '
                             'multi-level -> one-hot). Pass the sentinel 1000G_pops.txt (the default) '
                             'to skip. Port of the reference MakeCovariateFile.R behaviour.')

    args = parser.parse_args()

    print("=" * 60)
    print("Generating CovariatePCs.txt for sQTL analysis")
    print("=" * 60)

    # Read genotype PCs (keyed by genotype IDs)
    gen_pcs = read_genotype_pcs(args.genotype_pcs, args.num_gen_pcs)

    # Read splicing PCs. ProcessSplicing.R already renames PC sample IDs to genotype-ID space
    # (like eQTLGen ProcessExpression.R does for exp_PCs.txt), so no GTE remap is needed here.
    splicing_pcs = read_splicing_pcs(args.splicing_pcs, args.num_splicing_pcs, args.pc_prefix)

    # Merge genotype and splicing PCs (both keyed by genotype SampleID)
    print("\nMerging genotype and splicing PCs...")
    cov = pd.merge(gen_pcs, splicing_pcs, on='SampleID', how='inner')
    print(f"  After merge: {len(cov)} samples")

    # Check for sample mismatches
    gen_only = set(gen_pcs['SampleID']) - set(splicing_pcs['SampleID'])
    splicing_only = set(splicing_pcs['SampleID']) - set(gen_pcs['SampleID'])

    if gen_only:
        print(f"  WARNING: {len(gen_only)} samples only in genotype PCs")
    if splicing_only:
        print(f"  WARNING: {len(splicing_only)} samples only in splicing PCs (after GTE mapping)")

    # Add GenSex if provided
    if args.gen_sex:
        sex_df = read_gen_sex(args.gen_sex)
        if sex_df is not None:
            print("\nMerging GenSex...")
            cov = pd.merge(cov, sex_df, on='SampleID', how='left')
            print(f"  After merge: {len(cov)} samples")

            # Report missing GenSex
            missing_sex = cov['GenSex'].isna().sum()
            if missing_sex > 0:
                print(f"  WARNING: {missing_sex} samples have missing GenSex")
    else:
        print("\nNo GenSex file provided, skipping GenSex column")

    # Add cohort-specific additional covariates (port of MakeCovariateFile.R).
    # The sentinel 1000G_pops.txt means "none supplied" (matches the reference default).
    if args.additional_covariates and os.path.basename(args.additional_covariates) != '1000G_pops.txt':
        add_cov = read_additional_covariates(args.additional_covariates)
        missing = set(cov['SampleID']) - set(add_cov['SampleID'])
        if missing:
            sys.exit(
                f"ERROR: {len(missing)} sample(s) in the PC file are absent from the additional "
                f"covariate file. Every sample must be present in the covariate file. "
                f"Examples: {sorted(missing)[:5]}"
            )
        print("\nMerging additional covariates...")
        cov = pd.merge(cov, add_cov, on='SampleID', how='inner')
        print(f"  After merge: {len(cov)} samples, {cov.shape[1] - 1} covariate column(s)")
    else:
        print("\nNo additional covariates supplied (sentinel/none); skipping.")

    # Write output
    print(f"\nWriting output to: {args.output}")
    cov.to_csv(args.output, sep='\t', index=False, na_rep='NA')

    print(f"\nSummary:")
    print(f"  Total samples: {len(cov)}")
    print(f"  Genotype PCs: {args.num_gen_pcs}")
    print(f"  Splicing PCs: {args.num_splicing_pcs}")
    if 'GenSex' in cov.columns:
        print(f"  GenSex included: Yes")
    else:
        print(f"  GenSex included: No")

    print("\nColumn order:")
    print(f"  {' | '.join(cov.columns.tolist())}")

    print("\n" + "=" * 60)
    print("CovariatePCs.txt generation completed successfully!")
    print("=" * 60)

if __name__ == '__main__':
    main()
