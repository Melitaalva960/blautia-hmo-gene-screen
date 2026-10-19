## =====================================================================
## 04  Species-level heatmaps for the top N species (default 20)
##
##  A. Mean % identity per species x gene (absent genomes count as 0, so
##     a gene found in only a few genomes of a species appears pale)
##  B. Present/absent per species: a gene is "present" when at least
##     presence_frac (default 50%) of that species' genomes carry it
##  + a table of core / accessory / rare tiers (thesis definition)
##
## Input : results/filtering/identity_matrix.csv (from script 01), GTDB-Tk summary
## Output: results/figures/species_mean_identity_top20.png
##         results/figures/species_presence_50pct_top20.png
##         results/filtering/species_gene_summary.csv, gene_tiers.csv
## =====================================================================
source("R/utils.R")
suppressPackageStartupMessages({ library(pheatmap); library(RColorBrewer) })

s <- get_settings(list(
  matrix        = "results/filtering/identity_matrix.csv",
  taxonomy      = "results/taxonomy/gtdbtk.bac120.summary.tsv",
  gene_meta     = "data/panel/gene_metadata.tsv",
  out_dir       = "results/figures",
  table_dir     = "results/filtering",
  top_n         = 20,
  presence_frac = 0.5
))
ensure_dir(s$out_dir); ensure_dir(s$table_dir)

m   <- fread(s$matrix)
idm <- as.matrix(m[, -1, with = FALSE]); rownames(idm) <- m$genome
tax  <- read_taxonomy(s$taxonomy)
meta <- read_gene_metadata(s$gene_meta)

top <- top_species(tax, s$top_n)
tax_top <- tax[species %in% top & genome %in% rownames(idm)]
summ <- species_gene_summary(idm[tax_top$genome, , drop = FALSE], tax_top)
summ <- merge(summ, meta, by.x = "qseqid", by.y = "query", all.x = TRUE)
fwrite(summ, file.path(s$table_dir, "species_gene_summary.csv"))

tiers <- gene_tiers(summ, s$presence_frac, n_species = length(top))
tiers <- merge(tiers, meta, by.x = "qseqid", by.y = "query", all.x = TRUE)[order(-n_species_present)]
fwrite(tiers, file.path(s$table_dir, "gene_tiers.csv"))
print(tiers[, .N, by = tier])

## ---- Wide matrices (species x gene label) --------------------------------
to_wide <- function(value) {
  w <- dcast(summ, species ~ label, value.var = value, fill = 0)
  x <- as.matrix(w[, -1, with = FALSE]); rownames(x) <- short_species(w$species); x
}
mean_mat <- to_wide("mean_pident") / 100
pres_mat <- (to_wide("presence_frac") >= s$presence_frac) * 1

## Drop genes absent from every species and species with no variation
## (constant rows/columns break hierarchical clustering)
keep_cols <- apply(mean_mat, 2, sd) > 0
keep_rows <- apply(mean_mat[, keep_cols, drop = FALSE], 1, sd) > 0
mean_mat <- mean_mat[keep_rows, keep_cols, drop = FALSE]
pres_mat <- pres_mat[rownames(mean_mat), colnames(mean_mat), drop = FALSE]

## ---- Annotation colours ---------------------------------------------------
## Cluster the rows once, then colour species along the cluster order, so
## the species bar reads as a smooth gradient.
row_hc  <- hclust(dist(mean_mat), method = "ward.D2")
sp_cols <- setNames(colorRampPalette(c("#440154", "#3b528b", "#21908c", "#5dc963", "#fde725"))(nrow(mean_mat)),
                    rownames(mean_mat)[row_hc$order])
col_annot <- data.frame(Pathway = meta$pathway[match(colnames(mean_mat), meta$label)],
                        row.names = colnames(mean_mat))
pw <- sort(unique(col_annot$Pathway))
ann_colors <- list(Species = sp_cols,
                   Pathway = setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(pw)), pw))
row_annot <- data.frame(Species = rownames(mean_mat), row.names = rownames(mean_mat))

## Column dendrogram, mirrored left-to-right (same clusters, flipped order)
col_hc <- as.hclust(rev(as.dendrogram(hclust(dist(t(mean_mat)), method = "ward.D2"))))

common <- list(cluster_rows = row_hc, cluster_cols = col_hc,
               annotation_row = row_annot, annotation_col = col_annot,
               annotation_colors = ann_colors, fontsize_row = 10, fontsize_col = 7,
               angle_col = 90, width = 20, height = 10)

## ---- A. Mean % identity -----------------------------------------------------
f1 <- file.path(s$out_dir, sprintf("species_mean_identity_top%d.png", length(top)))
do.call(pheatmap, c(list(mean_mat, color = colorRampPalette(c("white", "#08306b"))(100),
                         breaks = seq(0, 1, length.out = 101),
                         main = sprintf("HMO gene presence across Blautia species (top %d by genome count): mean %% identity",
                                        length(top)),
                         filename = f1), common))

## ---- B. Present / absent (>= presence_frac of genomes) -------------------------
f2 <- file.path(s$out_dir, sprintf("species_presence_%.0fpct_top%d.png", 100 * s$presence_frac, length(top)))
do.call(pheatmap, c(list(pres_mat, color = c("#E7E9EC", "#0F6E56"), breaks = c(-0.5, 0.5, 1.5),
                         legend_breaks = c(0, 1), legend_labels = c("Absent", "Present"),
                         main = sprintf("HMO gene presence/absence (>= %.0f%% of a species' genomes): top %d species",
                                        100 * s$presence_frac, length(top)),
                         filename = f2), common))

for (f in c(f1, f2)) if (!file.exists(f) || file.info(f)$size == 0) stop("Empty output: ", f)
cat("Saved", f1, "and", f2, "\n")
