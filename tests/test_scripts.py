"""Unit tests for the filtering and matrix scripts.  Run with:  pytest -q"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from filter_hits import filter_hits, load_hits  # noqa: E402
from build_matrix import build_matrix, gene_tiers  # noqa: E402

RAW = (
    # qseqid sseqid pident length qlen qstart qend evalue bitscore
    "geneA\tG1|c1\t45\t300\t300\t1\t300\t1e-50\t400\n"   # pass
    "geneA\tG1|c2\t40\t300\t300\t1\t300\t1e-40\t300\n"   # pass, but lower bitscore
    "geneA\tG2|c1\t25\t300\t300\t1\t300\t1e-20\t150\n"   # fails identity (<30)
    "geneB\tG2|c1\t60\t150\t300\t1\t150\t1e-30\t200\n"   # fails coverage (50%)
    "geneB\tG3|c9\t35\t280\t300\t11\t290\t1e-3\t90\n"    # fails e-value
    "geneB\tG3|c1\t35\t280\t300\t11\t290\t1e-10\t120\n"  # pass (coverage 93.3%)
)


def write_raw(tmp_path):
    f = tmp_path / "raw.tsv"
    f.write_text(RAW)
    return f


def test_filter_keeps_best_hit_per_genome_gene(tmp_path):
    best = filter_hits(load_hits(write_raw(tmp_path)))
    pairs = set(zip(best.genome, best.qseqid))
    assert pairs == {("G1", "geneA"), ("G3", "geneB")}
    g1 = best[best.genome == "G1"].iloc[0]
    assert g1.sseqid == "G1|c1"          # highest bitscore kept
    assert g1.bitscore == 400


def test_coverage_calculation(tmp_path):
    hits = load_hits(write_raw(tmp_path))
    row = hits[hits.sseqid == "G3|c1"].iloc[0]
    assert abs(row.qcov - 93.333) < 0.01


def test_matrix_includes_genomes_and_genes_without_hits(tmp_path):
    best = filter_hits(load_hits(write_raw(tmp_path)))
    m = build_matrix(best, ["G1", "G2", "G3", "G4"], ["geneA", "geneB", "geneC"])
    assert m.shape == (4, 3)
    assert m.loc["G1", "geneA"] == 1
    assert m.loc["G2"].sum() == 0
    assert m["geneC"].sum() == 0


def test_tiers():
    m = pd.DataFrame({"core_gene": [1] * 20,
                      "acc_gene": [1] * 10 + [0] * 10,
                      "rare_gene": [1] + [0] * 19})
    t = gene_tiers(m).set_index("gene")["tier"]
    assert t["core_gene"] == "core"
    assert t["acc_gene"] == "accessory"
    assert t["rare_gene"] == "rare"
