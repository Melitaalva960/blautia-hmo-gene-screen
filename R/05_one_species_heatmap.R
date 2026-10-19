## =====================================================================
## 05  Genome-level heatmap for ONE species (default: B. wexlerae)
##
## Rows = every genome of that species, columns = genes coloured by
## pathway, cells = % identity of the best hit (white = absent).
## Shows variation in gene content within a species.
##
## Input : results/filtering/identity_matrix.csv (from script 01), GTDB-Tk summary
## Output: results/figures/heatmap_<species>_genomes.png
##
## Run:   Rscript R/05_one_species_heatmap.R species=wexlerae
## =====================================================================
source("R/utils.R")
suppressPackageStartupMessages({ library(pheatmap); library(RColorBrewer) })

s <- get_settings(list(
  matrix    = "results/filtering/identity_matrix.csv",
  taxonomy  = "results/taxonomy/gtdbtk.bac120.summary.tsv",
  gene_meta = "data/panel/gene_metadata.tsv",
  out_dir   = "results/figures",
  species   = "wexlerae"          # matched case-insensitively
))
ensure_dir(s$out_dir)

m   <- fread(s$matrix)
idm <- as.matrix(m[, -1, with = FALSE]); rownames(idm) <- m$genome
tax  <- read_taxonomy(s$taxonomy)
meta <- read_gene_metadata(s$gene_meta)

sel <- tax[str_detect(species, regex(s$species, ignore_case = TRUE)) & genome %in% rownames(idm)]
if (nrow(sel) == 0) stop("No genomes match species '", s$species, "'. Species found: ",
                         paste(head(unique(tax$species), 20), collapse = ", "))
cat(sprintf("%d genomes matched: %s\n", nrow(sel), paste(unique(sel$species), collapse = ", ")))

mat <- idm[sel$genome, , drop = FALSE] / 100
## Constant columns/rows (e.g. genes absent in every genome) break clustering
mat <- mat[, apply(mat, 2, sd) > 0, drop = FALSE]
mat <- mat[apply(mat, 1, sd) > 0, , drop = FALSE]
colnames(mat) <- meta$label[match(colnames(mat), meta$query)]

col_annot <- data.frame(Pathway = meta$pathway[match(colnames(mat), meta$label)],
                        row.names = colnames(mat))
pw <- sort(unique(col_annot$Pathway))
pw_cols <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(pw)), pw)

label <- short_species(unique(sel$species)[1])
out_file <- file.path(s$out_dir, sprintf("heatmap_%s_genomes.png", gsub("[^A-Za-z0-9]+", "_", s$species)))
pheatmap(mat,
         color = colorRampPalette(c("white", "#08306b"))(100), breaks = seq(0, 1, length.out = 101),
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = col_annot, annotation_colors = list(Pathway = pw_cols),
         show_rownames = FALSE, fontsize_col = 7, angle_col = 90,
         main = sprintf("HMO gene presence across all %s genomes (n = %d)", label, nrow(mat)),
         filename = out_file, width = 20, height = 14)
if (!file.exists(out_file) || file.info(out_file)$size == 0) stop("Empty output: ", out_file)
cat("Saved", out_file, "\n")
