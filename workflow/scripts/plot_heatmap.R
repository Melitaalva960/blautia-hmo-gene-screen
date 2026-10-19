#!/usr/bin/env Rscript
# Heatmaps of HMO / FOS / GOS gene presence across Blautia genomes.
#
# Usage:
#   Rscript plot_heatmap.R <matrix.tsv> <genome_heatmap.png> <species_heatmap.png> [top_n_species]
#
# 1. Genome-level heatmap: every genome (rows) x every panel gene (columns).
# 2. Species-level heatmap: % of genomes carrying each gene, for the top N species
#    by genome count (a gene is usually called "present" in a species at >= 50%).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: plot_heatmap.R matrix.tsv genome.png species.png [top_n]")
matrix_file <- args[1]
genome_png  <- args[2]
species_png <- args[3]
top_n       <- if (length(args) >= 4) as.integer(args[4]) else 20

mat <- read_tsv(matrix_file, show_col_types = FALSE)
if (!"species" %in% names(mat)) mat$species <- "unclassified"

long <- mat |>
  pivot_longer(-c(genome, species), names_to = "gene", values_to = "present")

# Order genes from most to least common so core genes sit on the left
gene_order <- long |>
  group_by(gene) |>
  summarise(prev = mean(present), .groups = "drop") |>
  arrange(desc(prev)) |>
  pull(gene)
long$gene <- factor(long$gene, levels = gene_order)

# ---- 1. Genome-level heatmap ----
genome_order <- long |>
  group_by(genome, species) |>
  summarise(n = sum(present), .groups = "drop") |>
  arrange(species, desc(n)) |>
  pull(genome)
long$genome <- factor(long$genome, levels = rev(genome_order))

p1 <- ggplot(long, aes(gene, genome, fill = factor(present))) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "grey92", "1" = "#1F5C7A"),
                    labels = c("Absent", "Present"), name = NULL) +
  labs(x = "Panel gene", y = "Genome",
       title = "Gene presence/absence across genomes") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        axis.text.y = element_blank(), panel.grid = element_blank())
ggsave(genome_png, p1, width = 12, height = 8, dpi = 300)

# ---- 2. Species-level heatmap (top N species) ----
top_species <- mat |> count(species, sort = TRUE) |> slice_head(n = top_n)

by_species <- long |>
  filter(species %in% top_species$species) |>
  group_by(species, gene) |>
  summarise(pct = 100 * mean(present), .groups = "drop") |>
  left_join(top_species, by = "species") |>
  mutate(label = sprintf("%s (n=%d)", species, n),
         label = reorder(label, n))

p2 <- ggplot(by_species, aes(gene, label, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "grey95", high = "#1F5C7A", limits = c(0, 100),
                      name = "% genomes") +
  labs(x = "Panel gene", y = NULL,
       title = sprintf("Gene prevalence in the top %d species", nrow(top_species))) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        panel.grid = element_blank())
ggsave(species_png, p2, width = 12, height = 6, dpi = 300)
