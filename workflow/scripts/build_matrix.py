#!/usr/bin/env python3
"""Turn filtered best hits into a binary gene presence/absence matrix.

Rows = genomes, columns = panel genes, 1 = present, 0 = absent.
Every genome in --genome-list appears in the matrix, even with no hits,
and every gene in --panel-genes appears, even if never found.

Optionally joins GTDB-Tk taxonomy (--taxonomy, a gtdbtk *.summary.tsv) and adds
a "species" column, then writes per-gene prevalence and a core/accessory/rare
classification:
    core      >= 95% of genomes
    accessory 15-95%
    rare      < 15%
"""
import argparse

import pandas as pd


def build_matrix(best_hits, genomes, genes):
    present = best_hits.assign(value=1).pivot_table(
        index="genome", columns="qseqid", values="value", aggfunc="max", fill_value=0)
    matrix = present.reindex(index=genomes, columns=genes, fill_value=0).astype(int)
    matrix.index.name = "genome"
    matrix.columns.name = None
    return matrix


def gene_tiers(matrix, core=95.0, rare=15.0):
    prevalence = matrix.mean(axis=0) * 100
    tier = pd.cut(prevalence, bins=[-0.1, rare - 1e-9, core - 1e-9, 100.1],
                  labels=["rare", "accessory", "core"])
    return pd.DataFrame({"gene": matrix.columns,
                         "prevalence_pct": prevalence.round(2).values,
                         "tier": tier.astype(str).values})


def species_from_gtdb(path):
    tax = pd.read_csv(path, sep="\t", usecols=["user_genome", "classification"])
    tax["species"] = tax["classification"].str.extract(r"s__([^;]*)$")[0].fillna("")
    tax.loc[tax["species"] == "", "species"] = "unclassified"
    return tax.set_index("user_genome")["species"]


def read_lines(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("best_hits", help="output of filter_hits.py")
    p.add_argument("--genome-list", required=True, help="one genome ID per line")
    p.add_argument("--panel-genes", required=True, help="one gene ID per line")
    p.add_argument("--taxonomy", help="GTDB-Tk summary TSV (optional)")
    p.add_argument("-o", "--output", required=True, help="presence/absence matrix TSV")
    p.add_argument("--tiers", required=True, help="per-gene prevalence and tier TSV")
    args = p.parse_args(argv)

    best = pd.read_csv(args.best_hits, sep="\t")
    matrix = build_matrix(best, read_lines(args.genome_list), read_lines(args.panel_genes))
    gene_tiers(matrix).to_csv(args.tiers, sep="\t", index=False)

    if args.taxonomy:
        species = species_from_gtdb(args.taxonomy)
        matrix.insert(0, "species", species.reindex(matrix.index).fillna("unclassified"))
    matrix.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
