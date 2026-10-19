## =====================================================================
## 01  Filter tBLASTn hits and build the presence/absence matrix
##
## Input : one tBLASTn file per genome (<genome>_tblastn.tsv, outfmt 6 with
##         qlen after qseqid; see R/utils.R)
## Output: results/filtering/
##           best_hits.csv               best hit per genome-gene pair
##           presence_absence_matrix.csv genomes x genes, 1/0
##           identity_matrix.csv         genomes x genes, % identity (0 = absent)
##           gene_prevalence.csv         % of genomes carrying each gene
##           summary.csv                 counts at each step
##
## Run:   Rscript R/01_filter_hits.R tblastn_dir=results/tblastn
## =====================================================================
source("R/utils.R")

s <- get_settings(list(
  tblastn_dir = "results/tblastn",
  gene_meta   = "data/panel/gene_metadata.tsv",
  out_dir     = "results/filtering",
  coverage    = "length",          # "length" (thesis) or "span"
  evalue_max  = THRESHOLDS$evalue_max,
  pident_min  = THRESHOLDS$pident_min,
  qcov_min    = THRESHOLDS$qcov_min
))
ensure_dir(s$out_dir)

## Every screened genome has a file, even if it had no hits, so the genome
## list comes from the file names (genomes without hits stay in the matrix).
files   <- list.files(s$tblastn_dir, pattern = "_tblastn\\.tsv$", full.names = TRUE)
genomes <- sort(genome_id_from_file(files))
cat("Genomes screened:", length(genomes), "\n")

hits <- read_tblastn_dir(s$tblastn_dir)
hits <- add_query_coverage(hits, s$coverage)
if (anyNA(hits$qcov)) warning("Some hits have no qlen, so coverage could not be computed; they are dropped.")

best <- filter_best_hits(hits, s$evalue_max, s$pident_min, s$qcov_min)

gene_meta <- read_gene_metadata(s$gene_meta)
unknown <- setdiff(unique(hits$qseqid), gene_meta$query)
if (length(unknown) > 0) warning("Query IDs not in gene metadata: ", paste(unknown, collapse = ", "))
genes <- union(gene_meta$query, unique(hits$qseqid))

pa  <- presence_matrix(best, genomes, genes)
idm <- identity_matrix(best, genomes, genes)

fwrite(best, file.path(s$out_dir, "best_hits.csv"))
fwrite(as.data.table(pa,  keep.rownames = "genome"), file.path(s$out_dir, "presence_absence_matrix.csv"))
fwrite(as.data.table(idm, keep.rownames = "genome"), file.path(s$out_dir, "identity_matrix.csv"))

prevalence <- data.table(query = colnames(pa),
                         n_genomes = colSums(pa),
                         pct_genomes = round(100 * colMeans(pa), 2))
prevalence <- merge(prevalence, gene_meta, by = "query", all.x = TRUE)[order(-pct_genomes)]
fwrite(prevalence, file.path(s$out_dir, "gene_prevalence.csv"))

summary_dt <- data.table(
  genomes_screened  = length(genomes),
  panel_genes       = length(genes),
  raw_hits          = nrow(hits),
  hits_passing      = nrow(hits[evalue <= s$evalue_max & pident >= s$pident_min & qcov >= s$qcov_min]),
  best_hits         = nrow(best),
  genes_detected    = uniqueN(best$qseqid),
  coverage_method   = s$coverage
)
fwrite(summary_dt, file.path(s$out_dir, "summary.csv"))
print(summary_dt)
cat("Saved results to", s$out_dir, "\n")
