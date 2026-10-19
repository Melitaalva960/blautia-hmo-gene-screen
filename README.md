# Blautia HMO / FOS / GOS gene screen

![tests](https://github.com/Melitaalva960/blautia-hmo-gene-screen/actions/workflows/tests.yml/badge.svg)

Comparative genomics pipeline from my MSc thesis at University College Cork, asking:
**can the gut bacterium *Blautia* use human milk oligosaccharides (HMOs) and the
FOS/GOS prebiotics added to infant formula?**

A custom 64-gene, 19-pathway HMO-utilisation panel was screened against **3,081 *Blautia*
genomes** with tBLASTn, and the results were summarised as gene presence/absence across
genomes and species.

## Why this project

*Blautia* is abundant in the human gut but understudied. Existing annotation tools are built
for *Bifidobacterium*. For example, **bifidoAnnotator annotated 18 glycoside hydrolases in
*B. breve* UCC2003 but 0 in *Blautia wexlerae* DSM 19850 and *Roseburia hominis***. This is a
limit of the tool, not proof that the genes are absent. A homology screen with a custom
panel detected **27 of the 64 panel genes in *B. wexlerae* DSM 19850**.

## Workflow

```
genome FASTA files (UHGG isolates + SPIRE MAGs)
   │
   ├─ CheckM2 ───────── keep completeness > 90%, contamination < 5%
   ├─ GTDB-Tk (r220) ── species assignment
   ├─ tBLASTn ───────── 64 panel proteins vs each genome (6-frame translation,
   │                    works on MAGs without gene calling)
   ├─ R/01 ──────────── E ≤ 1e-5, identity ≥ 30%, query coverage ≥ 70%,
   │                    best hit (bitscore) per genome-gene pair → matrices
   └─ R/02-07 ───────── species counts, heatmaps, core / accessory / rare tiers
```

## Key results (thesis)

| Gene set | Genomes carrying it |
|---|---|
| FOS-degradation enzymes (GH32) | **84.5%** |
| GOS-degradation enzymes (β-galactosidases) | **75.1%** |
| Dedicated GOS transporter (from *B. breve*) | **0.9%** |

- 74,359 filtered best hits covering 62 of the 64 HMO panel genes
- Genes fall into **core, accessory and rare tiers** across species, pointing to a division
  of labour and cross-feeding between *Blautia* species rather than complete HMO
  degradation by any single species
- Phylogeny of the HMO-utilisation gene RhGLnbp112 (MAFFT + IQ-TREE 2 with ModelFinder,
  rooted on *B. breve* UCC2003) placed *Roseburia hominis* inside the *Blautia* clade

<!-- Figures can be added to docs/ and shown here, e.g.
![Top 20 species](docs/species_presence_50pct_top20.png)
-->

## Analysis scripts (R)

Run from the repository root. Every setting can be changed on the command line as `key=value`.

| Script | What it does | Main output |
|---|---|---|
| `R/01_filter_hits.R` | Reads one tBLASTn file per genome, applies the thresholds, keeps the best hit per genome-gene pair | presence/absence and % identity matrices, gene prevalence |
| `R/02_species_counts.R` | Genomes per species from GTDB-Tk | bar charts (all species, top 20) |
| `R/03_genome_heatmap.R` | Every genome × gene, clustered, coloured by species and pathway | genome-level heatmap |
| `R/04_species_heatmaps.R` | Top 20 species: mean % identity, and present/absent (≥ 50% of a species' genomes); core/accessory/rare tiers | 2 heatmaps, tier table |
| `R/05_one_species_heatmap.R` | All genomes of one species (default *B. wexlerae*) | within-species heatmap |
| `R/06_pathway_group_heatmaps.R` | Fucose-based and LNT-family genes grouped by function, names coloured by tier | grouped heatmaps |
| `R/07_single_genome_identity.R` | % identity of all 64 genes in one genome | one-row heatmap |
| `R/utils.R` | Shared functions used by all scripts | |

**Tier definition:** a gene is present in a species if ≥ 50% of its genomes carry it. Across
the top 20 species, a gene is **core** in ≥ 18, **accessory** in 5–17 and **rare** in ≤ 4.

**Query coverage:** `coverage=length` (alignment length / query length) was used for the
thesis matrix; `coverage=span` ((qend − qstart + 1) / query length) is also available.
Both gave identical calls for *B. wexlerae* DSM 19850.

Example:

```bash
Rscript R/01_filter_hits.R tblastn_dir=results/tblastn
Rscript R/04_species_heatmaps.R taxonomy=results/taxonomy/gtdbtk.bac120.summary.tsv top_n=20
```

## Running the full workflow

```bash
conda env create -f environment.yml
conda activate blautia-screen

# genomes -> data/genomes/   panel proteins -> data/panel/hmo_panel.faa
# set the CheckM2 and GTDB-Tk database paths in config/config.yaml, then:
snakemake --cores 16 --use-conda

# quick run without CheckM2/GTDB-Tk (filtering only):
snakemake --cores 4 --use-conda --config run_qc=False
```

## Data

- `data/panel/gene_metadata.tsv`: the 64 panel genes with short labels and pathways
- `data/panel/gene_groups.tsv`: functional groups used for the Fucose and LNT figures
- Genomes are not stored here. They are public (UHGG, SPIRE); see `data/genomes/README.md`.

## Tests

```bash
Rscript tests/test_utils.R
```

Tests run automatically on GitHub Actions on every change.

## Tools

CheckM2 · GTDB-Tk · BLAST+ (tBLASTn, BLASTP) · MAFFT · IQ-TREE 2 · iTOL · Snakemake ·
R (data.table, dplyr, ggplot2, pheatmap, ggtext, patchwork) · CAZy · UHGG · SPIRE

## Limitations

- Homology at ≥ 30% identity suggests, but does not prove, enzyme function. Growth
  experiments are needed to confirm utilisation.
- Genes split across fragmented MAG contigs can fall below the coverage threshold, so
  some genes may be missed.

## Author

**Alva Melita Alwin** · MSc Bioinformatics & Computational Biology, University College Cork  
Supervisor: Dr Hilary Browne, School of Microbiology, UCC  
Open to bioinformatics and genomics roles in Ireland
