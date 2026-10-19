## =====================================================================
## 03  Genome-level presence/absence heatmap (all genomes x panel genes)
##
## Rows = genomes (clustered, coloured by species), columns = genes
## (clustered, coloured by pathway). Species outside the top N are
## grouped as "Other" so the colour bar stays readable.
##
## Input : results/filtering/presence_absence_matrix.csv (from script 01)
##         GTDB-Tk summary
## Output: results/figures/heatmap_all_genomes.png
## =====================================================================
source("R/utils.R")
suppressPackageStartupMessages({ library(pheatmap); library(RColorBrewer) })

s <- get_settings(list(
  matrix    = "results/filtering/presence_absence_matrix.csv",
  taxonomy  = "results/taxonomy/gtdbtk.bac120.summary.tsv",
  gene_meta = "data/panel/gene_metadata.tsv",
  out_dir   = "results/figures",
  top_n     = 20
))
ensure_dir(s$out_dir)

m   <- fread(s$matrix)
mat <- as.matrix(m[, -1, with = FALSE]); rownames(mat) <- m$genome
tax <- read_taxonomy(s$taxonomy)
meta <- read_gene_metadata(s$gene_meta)

## Genes never detected carry no information for clustering
mat <- mat[, colSums(mat) > 0, drop = FALSE]
colnames(mat) <- meta$label[match(colnames(mat), meta$query)]

## Row annotation: species (top N, else "Other")
keep <- top_species(tax, s$top_n)
sp <- tax$species[match(rownames(mat), tax$genome)]
sp[is.na(sp)] <- "Unclassified"
sp <- ifelse(sp %in% keep, short_species(sp), "Other")
row_annot <- data.frame(Species = sp, row.names = rownames(mat))

col_annot <- data.frame(Pathway = meta$pathway[match(colnames(mat), meta$label)],
                        row.names = colnames(mat))

sp_levels <- sort(setdiff(unique(sp), "Other"))
sp_cols <- setNames(colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(length(sp_levels)),
                    sp_levels)
if ("Other" %in% sp) sp_cols <- c(sp_cols, Other = "grey80")
pw <- sort(unique(col_annot$Pathway))
pw_cols <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(pw)), pw)

out_file <- file.path(s$out_dir, "heatmap_all_genomes.png")
pheatmap(mat,
         color = c("white", "#08306b"), breaks = c(-0.5, 0.5, 1.5),
         legend_breaks = c(0, 1), legend_labels = c("Absent", "Present"),
         clustering_distance_rows = "binary", clustering_distance_cols = "binary",
         clustering_method = "complete",
         annotation_row = row_annot, annotation_col = col_annot,
         annotation_colors = list(Species = sp_cols, Pathway = pw_cols),
         show_rownames = FALSE, fontsize_col = 7, angle_col = 90,
         main = sprintf("HMO gene presence/absence in %d Blautia genomes", nrow(mat)),
         filename = out_file, width = 20, height = 14)
cat("Saved", out_file, "\n")
