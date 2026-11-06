#!/usr/bin/env python3
"""
transpose_matrix_np.py

Memory-efficient transpose for a WIDE matrix (few rows, MANY columns), e.g. the samples x ~670k-feature
splicing phenotype. `pandas.read_csv` on hundreds of thousands of COLUMNS has enormous per-column overhead
and OOMs (the raw data is only ~24 GB); numpy reads the same data into a compact 2-D array with no
per-column overhead. The transpose is then TALL (many rows, few columns), which pandas writes quickly and
safely. Drop-in for transpose_matrix.py's interface (--input/--output/--junction_list).

Input : row id (col 1) + value columns; first line = header (id-name + column names).
Output: the transpose (input columns become output rows). --junction_list gets the input ROW ids
        (= output column headers), mirroring transpose_matrix.py.
"""

import argparse
import sys
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--junction_list', required=True)
    args = ap.parse_args()

    with open(args.input) as f:
        header = f.readline().rstrip('\n').split('\t')
    col_names = header[1:]                      # input columns -> output rows
    n_cols = len(col_names)

    with open(args.input) as f:                 # count data rows (cheap pass)
        f.readline()
        n_rows = sum(1 for _ in f)
    print(f"Reading {n_rows} rows x {n_cols} cols with numpy...", file=sys.stderr)

    A = np.empty((n_rows, n_cols), dtype=np.float64)
    row_ids = []
    with open(args.input) as f:
        f.readline()
        for i, line in enumerate(f):
            p = line.rstrip('\n').split('\t')
            row_ids.append(p[0])
            A[i, :] = np.asarray(p[1:], dtype=np.float64)

    with open(args.junction_list, 'w') as g:    # input row ids = output column headers
        g.write('\n'.join(row_ids) + '\n')

    # Output is TALL (n_cols rows x n_rows cols) -> pandas handles the few columns efficiently.
    print("Writing transpose with pandas...", file=sys.stderr)
    df = pd.DataFrame(A.T, index=col_names, columns=row_ids)
    df.index.name = 'ID'
    df.to_csv(args.output, sep='\t')
    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
