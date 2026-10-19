# Blautia HMO / FOS / GOS gene screen

![tests](https://github.com/Melitaalva960/blautia-hmo-gene-screen/actions/workflows/tests.yml/badge.svg)

A reproducible **Snakemake** pipeline for screening thousands of bacterial genomes
for carbohydrate-utilisation genes. It was built for my MSc thesis at University College Cork,
which asked: **can the gut bacterium *Blautia* use human milk oligosaccharides (HMOs)
and the FOS/GOS prebiotics added to infant formula?**

## Why this project

*Blautia* is abundant in the human gut but understudied. Existing annotation tools
were built for *Bifidobacterium*. For example, **bifidoAnnotator annotated 18 glycoside
hydrolases in *B. breve* UCC2003 but 0 in *Blautia wexlerae* DSM 19850 and *Roseburia hominis***.
This reflects the tool's limits, not a true absence of the genes. A homology-based screen
with a custom gene panel solves this.

## Workflow

```
genome FASTA files
   │
   ├─ CheckM2 ──────── keep completeness > 90%, contamination < 5%
   ├─ GTDB-Tk (r220) ─ species assignment
   │
   ├─ tBLASTn ──────── protein gene panel vs all genomes (6-frame translation,
   │                   works on MAGs without gene calling)
   ├─ filter ───────── E ≤ 1e-5, identity ≥ 30%, query coverage ≥ 70%,
   │                   best hit (bitscore) per genome-gene pair
   ├─ matrix ───────── binary presence/absence + core / accessory / rare tiers
   └─ R / ggplot2 ──── genome-level and species-level heatmaps
```

## Key results (MSc thesis, 3,081 *Blautia* genomes)

| Gene set | Genomes carrying it |
|---|---|
| FOS-degradation enzymes (GH32) | **84.5%** |
| GOS-degradation enzymes (β-galactosidases) | **75.1%** |
| Dedicated GOS transporter (from *B. breve*) | **0.9%** |

- 74,359 filtered best hits covering 62 of the 64 HMO panel genes
- HMO genes fall into **core, accessory and rare tiers**, pointing to a division of labour
  and cross-feeding between *Blautia* species rather than complete degradation by one species
- Phylogeny of the HMO-utilisation gene RhGLnbp112 (MAFFT + IQ-TREE 2 with ModelFinder)
  placed *Roseburia hominis* inside the *Blautia* clade, with *B. breve* UCC2003 as the outgroup

<!-- Add figures here once approved, e.g.
![Species heatmap](docs/heatmap_species.png)
-->

## Quick start

```bash
# 1. Get the code
git clone https://github.com/Melitaalva960/blautia-hmo-gene-screen.git
cd blautia-hmo-gene-screen

# 2. Create the environment
conda env create -f environment.yml
conda activate blautia-screen

# 3. Add data
#    - genome FASTA files  -> data/genomes/   (see data/genomes/README.md)
#    - gene panel proteins -> data/panel/hmo_panel.faa

# 4. Quick test (skips CheckM2/GTDB-Tk, which need large databases)
snakemake --cores 4 --use-conda --config run_qc=False

# 5. Full run (set database paths in config/config.yaml first)
snakemake --cores 16 --use-conda
```

Outputs appear in `results/`:

| File | Content |
|---|---|
| `results/blast/best_hits.tsv` | filtered best hit per genome-gene pair |
| `results/matrix/presence_absence.tsv` | genomes × genes (1/0), with species |
| `results/matrix/gene_tiers.tsv` | prevalence and core/accessory/rare tier per gene |
| `results/figures/heatmap_*.png` | heatmaps |

All thresholds are in `config/config.yaml`.

## Tests

```bash
pytest -q
```

Tests run automatically on every push with GitHub Actions.

## Tools

Snakemake · CheckM2 · GTDB-Tk · BLAST+ (tBLASTn, BLASTP) · MAFFT · IQ-TREE 2 · iTOL ·
Python (pandas) · R (tidyverse, ggplot2) · CAZy · UHGG · SPIRE

## Limitations

- Homology at ≥ 30% identity suggests, but does not prove, enzyme function; growth
  experiments are needed to confirm utilisation.
- Genes on fragmented MAG contigs can fall below the coverage threshold, so some
  genes may be missed.

## Author

**Alva Melita Alwin** · MSc Bioinformatics & Computational Biology, University College Cork  
Supervisor: Dr Hilary Browne, School of Microbiology, UCC  
[LinkedIn](https://www.linkedin.com/in/YOUR-LINKEDIN) · Open to bioinformatics and genomics roles in Ireland
