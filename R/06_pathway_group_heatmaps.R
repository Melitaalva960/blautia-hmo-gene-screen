## =====================================================================
## 06  Pathway-group heatmaps (e.g. Fucose-based genes, LNT-family genes)
##
## For each group in data/panel/gene_groups.tsv: top N species (columns)
## x genes (rows, grouped by function), green = present in >= 50% of the
## species' genomes. Gene names are coloured by prevalence tier, computed
## from the data (green = core, amber = accessory, red = rare).
##
## Input : results/filtering/species_gene_summary.csv (from script 04)
## Output: results/figures/group_<group>_heatmap.png / .pdf
## =====================================================================
source("R/utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(ggtext) })

s <- get_settings(list(
  summary       = "results/filtering/species_gene_summary.csv",
  groups        = "data/panel/gene_groups.tsv",
  out_dir       = "results/figures",
  presence_frac = 0.5
))
ensure_dir(s$out_dir)

summ   <- fread(s$summary)
groups <- fread(s$groups, sep = "\t", quote = "")
missing <- setdiff(groups$label, summ$label)
if (length(missing) > 0) warning("Genes in gene_groups.tsv not found in the summary: ",
                                 paste(missing, collapse = ", "))

n_sp <- uniqueN(summ$species)
species_order <- summ[, .(n = max(n_genomes)), by = species][order(-n), species]
tiers <- gene_tiers(summ, s$presence_frac, n_sp)
tiers <- merge(tiers, unique(summ[, .(qseqid, label)]), by = "qseqid")
tier_colors <- c(core = "#0F6E56", accessory = "#B8860B", rare = "#99231D")

for (g in unique(groups$group)) {
  gg <- groups[group == g]
  cats <- unique(gg$category)                 # category order as written in the file
  gg[, category := factor(category, levels = cats)]
  setorder(gg, category)                      # keeps file order within each category

  d <- summ[label %in% gg$label, .(species, label, present = presence_frac >= s$presence_frac)]
  d <- merge(d, gg[, .(label, category)], by = "label")
  d[, label := factor(label, levels = rev(gg$label))]
  d[, species := factor(short_species(species), levels = short_species(species_order))]

  tier_of <- setNames(tiers$tier, tiers$label)
  y_labels <- setNames(sprintf("<span style='color:%s'>%s</span>",
                               tier_colors[tier_of[levels(d$label)]], levels(d$label)),
                       levels(d$label))

  p <- ggplot(d, aes(species, label, fill = present)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    facet_grid(rows = vars(category), scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c(`TRUE` = "#0F6E56", `FALSE` = "#E3E1D9"),
                      labels = c(`TRUE` = "Present", `FALSE` = "Absent"), name = NULL) +
    scale_y_discrete(labels = y_labels) +
    scale_x_discrete(position = "top") +
    labs(title = sprintf("%s HMO genes (%d)", g, nrow(gg)),
         subtitle = sprintf("Present in >= %.0f%% of a species' genomes, top %d Blautia species. Grouped by function: %s",
                            100 * s$presence_frac, n_sp, paste(cats, collapse = " → ")),
         caption = sprintf("Gene names coloured by prevalence tier: green = core (>= %d/%d species), amber = accessory (%d-%d), red = rare (<= %d).",
                           ceiling(0.9 * n_sp), n_sp, ceiling(0.25 * n_sp), ceiling(0.9 * n_sp) - 1,
                           ceiling(0.25 * n_sp) - 1)) +
    theme_minimal(base_family = "serif") +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          plot.caption = element_text(size = 7.5, colour = "grey45", hjust = 0),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 8, face = "bold", angle = 45, hjust = 0),
          axis.text.y = element_markdown(size = 9, face = "bold", hjust = 1),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(face = "bold", angle = 0, size = 9),
          strip.background = element_blank(),
          legend.position = "bottom",
          panel.spacing.y = unit(0.4, "lines"))

  base <- file.path(s$out_dir, sprintf("group_%s_heatmap", gsub("[^A-Za-z0-9]+", "_", g)))
  h <- 2.5 + 0.25 * nrow(gg)
  ggsave(paste0(base, ".png"), p, width = 11, height = h, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), p, width = 11, height = h)
  cat("Saved", paste0(base, ".png"), "\n")
}
