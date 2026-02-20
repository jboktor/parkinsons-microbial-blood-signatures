#!/usr/bin/env python3
"""
Run clusTCR clustering from the command line.

Example:
    python run_clusTCR.py \
        --tcrb /path/to/amppd_tcrb_2025-08-04.tsv \
        --outdir /path/to/results \
        --outname clusTCR_vgene_restricted \
        --include-vgene \
        --cdr3-col CDR3b_aa \
        --vgene-col TRBV \
        --n-cpus all
"""

import argparse
import sys
import pickle
from pathlib import Path
from datetime import datetime
from typing import Union, Optional

import pandas as pd
from clustcr import Clustering


def infer_sep(path: Path, explicit_sep: Optional[str]) -> Optional[str]:
    if explicit_sep is not None:
        return explicit_sep
    suf = path.suffix.lower()
    if suf in {".tsv", ".txt"}:
        return "\t"
    if suf in {".csv"}:
        return ","
    # Let pandas guess (e.g., for parquet we'll use a different reader)
    return None


def load_table(path: Path, sep: Optional[str]) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    return pd.read_table(path, sep=sep) if sep else pd.read_table(path)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run clusTCR clustering on a TCRB table."
    )
    p.add_argument("--tcrb", required=True, type=Path,
                   help="Path to input TCRB table (.tsv/.csv/.parquet).")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Directory to write outputs into (created if missing).")
    p.add_argument("--outname", required=True, type=str,
                   help="Base output file name (no extension).")
    p.add_argument("--include-vgene", action="store_true",
                   help="Include V-gene restriction during clustering.")
    p.add_argument("--cdr3-col", default="CDR3b_aa",
                   help="Column name for CDR3 amino-acid sequence (default: CDR3b_aa).")
    p.add_argument("--vgene-col", default="TRBV",
                   help="Column name for V gene (default: TRBV). Used only if --include-vgene.")
    p.add_argument("--n-cpus", default="all",
                   help="Number of CPUs to use (integer) or 'all' (default).")
    p.add_argument("--sep", default=None,
                   help="Field delimiter for CSV/TSV (overrides auto-detect).")
    p.add_argument("--no-pickle", action="store_true",
                   help="Do not write the pickle output.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # Validate paths
    if not args.tcrb.exists():
        print(f"[ERROR] Input file not found: {args.tcrb}", file=sys.stderr)
        return 2
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Determine separator & load
    sep = infer_sep(args.tcrb, args.sep)
    print(f"[INFO] Loading table from: {args.tcrb}")
    try:
        cdr3_df = load_table(args.tcrb, sep)
    except Exception as e:
        print(f"[ERROR] Failed to read input: {e}", file=sys.stderr)
        return 3

    # Check required columns
    if args.cdr3_col not in cdr3_df.columns:
        print(f"[ERROR] Missing CDR3 column '{args.cdr3_col}' in input.", file=sys.stderr)
        return 4
    if args.include_vgene and args.vgene_col not in cdr3_df.columns:
        print(f"[ERROR] --include-vgene set, but V-gene column '{args.vgene_col}' is missing.", file=sys.stderr)
        return 5

    # Instantiate Clustering
    n_cpus = args.n_cpus
    if isinstance(n_cpus, str) and n_cpus.lower() != "all":
        try:
            n_cpus = int(n_cpus)
        except ValueError:
            print("[ERROR] --n-cpus must be an integer or 'all'.", file=sys.stderr)
            return 6

    print(f"[INFO] Initializing Clustering(n_cpus={n_cpus})")
    clustering = Clustering(n_cpus=n_cpus)

    # Fit
    print(f"[INFO] Running clustering (include_vgene={args.include_vgene})...")
    try:
        output = clustering.fit(
            cdr3_df,
            include_vgene=args.include_vgene,
            cdr3_col=args.cdr3_col,
            v_gene_col=args.vgene_col if args.include_vgene else None,
        )
    except Exception as e:
        print(f"[ERROR] clusTCR.fit failed: {e}", file=sys.stderr)
        return 7

    # Prepare output paths
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base = f"{args.outname}"
    clusters_csv = args.outdir / f"{base}_clusters.csv"
    summary_csv = args.outdir / f"{base}_summary.csv"
    pkl_path = args.outdir / f"{base}.pkl"

    # Write outputs
    try:
        print(f"[INFO] Writing clusters to: {clusters_csv}")
        output.clusters_df.to_csv(clusters_csv, index=True)

        print(f"[INFO] Writing summary to: {summary_csv}")
        output.summary().to_csv(summary_csv, index=True)

        if not args.no_pickle:
            print(f"[INFO] Writing pickle to: {pkl_path}")
            with open(pkl_path, "wb") as f:
                pickle.dump(output, f)
    except Exception as e:
        print(f"[ERROR] Failed to write outputs: {e}", file=sys.stderr)
        return 8

    print("[INFO] Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())