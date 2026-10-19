## =====================================================================
## 02  Genome counts per species (all species, and top 20)
##
## Input : GTDB-Tk summary (user_genome, classification)
## Output: results/figures/species_counts_all.png
##         results/figures/species_counts_top20.png
##
## Run:   Rscript R/02_species_counts.R taxonomy=results/taxonomy/gtdbtk.bac120.summary.tsv
## =====================================================================
source("R/utils.R")
suppressPackageStartupMessages(library(ggplot2))

s <- get_settings(list(
  taxonomy = "results/taxonomy/gtdbtk.bac120.summary.tsv",
  out_dir  = "results/figures",
  top_n    = 20
))
ensure_dir(s$out_dir)

tax <- read_taxonomy(s$taxonomy)
counts <- tax[, .(n_genomes = .N), by = .(genus, species)][order(-n_genomes)]
cat(sprintf("%d genomes, %d genera, %d species\n",
            nrow(tax), uniqueN(tax$genus), uniqueN(tax$species)))

count_plot <- function(d, title, subtitle, text_size, base) {
  d <- copy(d)[, species := factor(species, levels = rev(species))]
  ggplot(d, aes(species, n_genomes, fill = genus)) +
    geom_col() +
    geom_text(aes(label = n_genomes), hjust = -0.15, size = text_size, fontface = "bold") +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Number of genomes", fill = "Genus") +
    theme_minimal(base_size = base) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"),
          axis.text.y = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "top", legend.justification = "left",
          legend.title = element_text(face = "bold"),
          plot.margin = margin(15, 55, 15, 15))
}

p_all <- count_plot(counts, "Blautia genome collection: genomes per species",
                    sprintf("%d genomes total", sum(counts$n_genomes)), 3.2, 12)
ggsave(file.path(s$out_dir, "species_counts_all.png"), p_all, width = 11,
       height = 0.22 * nrow(counts) + 2, limitsize = FALSE, dpi = 200, bg = "white")

top <- counts[seq_len(min(s$top_n, .N))]
p_top <- count_plot(top, sprintf("Top %d Blautia species by genome count", nrow(top)),
                    sprintf("%d of %d genomes", sum(top$n_genomes), sum(counts$n_genomes)), 4, 13)
ggsave(file.path(s$out_dir, sprintf("species_counts_top%d.png", nrow(top))), p_top,
       width = 13, height = 9.5, dpi = 250, bg = "white")

cat("Saved figures to", s$out_dir, "\n")
