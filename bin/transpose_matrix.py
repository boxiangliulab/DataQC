#!/usr/bin/env python3
"""
Transpose splicing matrix from junctions × samples to samples × junctions
to match eQTLGen expression format
"""

import argparse
import pandas as pd
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Transpose matrix to eQTLGen format')
    parser.add_argument('--input', required=True, help='Input matrix (junctions × samples)')
    parser.add_argument('--output', required=True, help='Output matrix (samples × junctions)')
    parser.add_argument('--junction_list', required=True, help='Output junction list file')
    return parser.parse_args()

def main():
    args = parse_args()

    print(f"Reading matrix from {args.input}...", file=sys.stderr)
    # Input format: junction | sample1 | sample2 | ...
    df = pd.read_csv(args.input, sep='\t', index_col=0)

    n_junctions, n_samples = df.shape
    print(f"Input: {n_junctions} junctions × {n_samples} samples", file=sys.stderr)

    # Save junction list (row names from input = column names in output)
    junction_ids = df.index.tolist()
    with open(args.junction_list, 'w') as f:
        for junc in junction_ids:
            f.write(f"{junc}\n")
    print(f"Saved {len(junction_ids)} junctions to {args.junction_list}", file=sys.stderr)

    # Transpose: samples × junctions
    print("Transposing matrix...", file=sys.stderr)
    df_transposed = df.T

    # Reset index to make sample IDs a column named "ID" (eQTLGen format)
    df_transposed.index.name = 'ID'

    print(f"Output: {df_transposed.shape[0]} samples × {df_transposed.shape[1]} junctions", file=sys.stderr)

    # Save transposed matrix
    print(f"Writing transposed matrix to {args.output}...", file=sys.stderr)
    df_transposed.to_csv(args.output, sep='\t', index=True)

    print("Done!", file=sys.stderr)

if __name__ == '__main__':
    main()
