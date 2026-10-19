## =====================================================================
## 07  % identity of every panel gene in ONE genome
##     (thesis example: Blautia wexlerae DSM 19850)
##
## A one-row heatmap of best-hit % identity for all 64 panel genes, with
## a coloured strip showing each gene's pathway. Useful for comparing a
## single reference genome against the panel.
##
## Input : one tBLASTn file (12- or 13-column outfmt 6) and the panel FASTA
##         (used for query lengths when the file has no qlen column)
## Output: results/figures/single_genome_identity.png
##
## Run:   Rscript R/07_single_genome_identity.R tblastn=path/B_wexlerae_tblastn.tsv \
##            panel=data/panel/hmo_panel.faa name="B. wexlerae DSM 19850"
## =====================================================================
source("R/utils.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

s <- get_settings(list(
  tblastn   = "",
  panel     = "data/panel/hmo_panel.faa",
  gene_meta = "data/panel/gene_metadata.tsv",
  name      = "B. wexlerae",
  out_file  = "results/figures/single_genome_identity.png",
  coverage  = "span"
))
if (s$tblastn == "") stop("Give the tBLASTn file: tblastn=path/to/file.tsv")
ensure_dir(dirname(s$out_file))

## Query lengths from the panel FASTA (needed for 12-column files)
fasta_lengths <- function(path) {
  lines <- readLines(path, warn = FALSE)
  h <- which(startsWith(lines, ">"))
  ends <- c(h[-1] - 1, length(lines))
  setNames(mapply(function(a, b) if (a > b) 0L else sum(nchar(lines[a:b])), h + 1, ends),
           trimws(sub("^>", "", lines[h])))
}

hits <- read_tblastn_file(s$tblastn)
if (is.null(hits)) stop("No hits read from ", s$tblastn)
if (all(is.na(hits$qlen))) {
  if (!file.exists(s$panel)) stop("File has no qlen column; give the panel FASTA with panel=...")
  hits[, qlen := as.numeric(fasta_lengths(s$panel)[qseqid])]
}
hits <- add_query_coverage(hits, s$coverage)
best <- filter_best_hits(hits)

meta <- read_gene_metadata(s$gene_meta)
pathway_order <- unique(meta$pathway)
d <- merge(meta, best[, .(query = qseqid, pident)], by = "query", all.x = TRUE)
d[is.na(pident), pident := 0]
d[, pathway := factor(pathway, levels = pathway_order)]
setorder(d, pathway, label)
d[, idx := .I]

palette <- setNames(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
                      "#e5c494", "#b3b3b3", "#1f78b4", "#33a02c", "#fb9a99", "#fdbf6f",
                      "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#a6cee3", "#f781bf",
                      "#999999")[seq_along(pathway_order)], pathway_order)

p_bar <- ggplot(d, aes(idx, 1, fill = pathway)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = palette, name = "Pathway", drop = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))

p_heat <- ggplot(d, aes(idx, s$name, fill = pident)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#08306b", limits = c(0, 100), name = "% identity\n(best hit)") +
  scale_x_continuous(breaks = d$idx, labels = d$label, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(face = "bold", size = 11),
        panel.grid = element_blank())

final <- (p_bar / p_heat) +
  plot_layout(heights = c(1, 6), guides = "collect") +
  plot_annotation(
    title = sprintf("HMO gene presence in %s (tBLASTn %% identity)", s$name),
    subtitle = sprintf("Filters: e-value <= %.0e, %% identity >= %d, query coverage >= %d%%. %d of %d genes detected.",
                       THRESHOLDS$evalue_max, THRESHOLDS$pident_min, THRESHOLDS$qcov_min,
                       sum(d$pident > 0), nrow(d)),
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey30")))

ggsave(s$out_file, final, width = 20, height = 7.5, dpi = 200, bg = "white")
cat("Saved", s$out_file, "\n")
