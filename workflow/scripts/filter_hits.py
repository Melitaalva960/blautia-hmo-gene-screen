#!/usr/bin/env python3
"""Filter tBLASTn/BLASTP hits and keep the best hit per genome-gene pair.

Input is BLAST tabular output (-outfmt 6) with these columns, in this order:
    qseqid stitle pident length qlen qstart qend evalue bitscore

The subject title (2nd column) must start with "GENOMEID|contig_name" so each hit
can be traced back to its genome (the Snakefile adds this prefix when building
the database). It is stored in the "sseqid" column below.

Filters (defaults match the thesis):
    E-value <= 1e-5, percent identity >= 30, query coverage >= 70%
Then, for each genome-gene pair, only the hit with the highest bitscore is kept.
"""
import argparse
import sys

import pandas as pd

COLUMNS = ["qseqid", "sseqid", "pident", "length", "qlen",
           "qstart", "qend", "evalue", "bitscore"]


def load_hits(path):
    hits = pd.read_csv(path, sep="\t", header=None, names=COLUMNS)
    if hits.empty:
        return hits.assign(genome=pd.Series(dtype=str), qcov=pd.Series(dtype=float))
    hits["genome"] = hits["sseqid"].str.split("|", n=1).str[0]
    # Query coverage = aligned part of the query / query length
    hits["qcov"] = (hits["qend"] - hits["qstart"] + 1) / hits["qlen"] * 100
    return hits


def filter_hits(hits, max_evalue=1e-5, min_identity=30.0, min_coverage=70.0):
    kept = hits[
        (hits["evalue"] <= max_evalue)
        & (hits["pident"] >= min_identity)
        & (hits["qcov"] >= min_coverage)
    ]
    # Best hit (highest bitscore) per genome-gene pair
    best = (kept.sort_values("bitscore", ascending=False)
                .drop_duplicates(subset=["genome", "qseqid"], keep="first")
                .sort_values(["genome", "qseqid"])
                .reset_index(drop=True))
    return best


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("hits", help="BLAST tabular output (outfmt 6, 9 columns)")
    p.add_argument("-o", "--output", required=True, help="filtered best hits (TSV)")
    p.add_argument("--max-evalue", type=float, default=1e-5)
    p.add_argument("--min-identity", type=float, default=30.0)
    p.add_argument("--min-coverage", type=float, default=70.0)
    args = p.parse_args(argv)

    hits = load_hits(args.hits)
    best = filter_hits(hits, args.max_evalue, args.min_identity, args.min_coverage)
    best.to_csv(args.output, sep="\t", index=False)
    print(f"{len(hits)} raw hits -> {len(best)} filtered best hits", file=sys.stderr)


if __name__ == "__main__":
    main()
