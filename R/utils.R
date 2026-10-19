## =====================================================================
## Shared helper functions used by all analysis scripts in R/
##
## Every script starts with:   source("R/utils.R")
## and is run from the repository root, e.g.
##   Rscript R/01_filter_hits.R tblastn_dir=results/tblastn
## =====================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

## ---- Command-line settings ------------------------------------------
## Each script defines its defaults; any of them can be overridden on the
## command line as key=value, e.g.  Rscript R/04_species_heatmaps.R top_n=10
get_settings <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  for (a in args) {
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2) stop("Arguments must look like key=value, got: ", a)
    key <- kv[1]; value <- paste(kv[-1], collapse = "=")
    if (!key %in% names(defaults)) {
      stop("Unknown setting '", key, "'. Valid settings: ",
           paste(names(defaults), collapse = ", "))
    }
    defaults[[key]] <- if (is.numeric(defaults[[key]])) as.numeric(value) else value
  }
  defaults
}

## ---- Thesis filter thresholds ---------------------------------------
THRESHOLDS <- list(evalue_max = 1e-5, pident_min = 30, qcov_min = 70)

## ---- Reading tblastn output -----------------------------------------
## The thesis tblastn files have 13 columns (outfmt 6 with qlen added
## after qseqid). Standard 12-column outfmt 6 is also accepted, but then
## there is no qlen and query coverage cannot be computed.
BLAST_COLS_13 <- c("qseqid", "qlen", "sseqid", "pident", "length", "mismatch",
                   "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
BLAST_COLS_12 <- setdiff(BLAST_COLS_13, "qlen")
NUMERIC_COLS  <- c("qlen", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore")

## "GUT_GENOME000123_tblastn.tsv" -> "GUT_GENOME000123"
genome_id_from_file <- function(path) {
  str_remove(tools::file_path_sans_ext(basename(path)), "_tblastn.*$|_blast.*$")
}

read_tblastn_file <- function(path) {
  if (file.size(path) == 0) return(NULL)            # genome with no hits
  dt <- tryCatch(fread(path, sep = "\t", header = FALSE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  if (ncol(dt) == 13) {
    setnames(dt, BLAST_COLS_13)
  } else if (ncol(dt) == 12) {
    setnames(dt, BLAST_COLS_12)
    dt[, qlen := NA_real_]
  } else {
    warning(sprintf("Skipping %s: %d columns (expected 12 or 13).", basename(path), ncol(dt)))
    return(NULL)
  }
  ## Force numeric types: one odd file would otherwise turn a whole
  ## combined column into character.
  dt[, (NUMERIC_COLS) := lapply(.SD, as.numeric), .SDcols = NUMERIC_COLS]
  dt[, genome := genome_id_from_file(path)]
  dt[]
}

## Read every *_tblastn.tsv file in a folder (optionally only some genomes)
read_tblastn_dir <- function(dir, genomes = NULL, pattern = "_tblastn\\.tsv$") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!is.null(genomes)) files <- files[genome_id_from_file(files) %in% genomes]
  if (length(files) == 0) stop("No files matching '", pattern, "' in ", dir)
  dt <- rbindlist(lapply(files, read_tblastn_file), use.names = TRUE)
  if (nrow(dt) == 0) stop("No hits could be read from ", dir)
  dt
}

## ---- Query coverage ---------------------------------------------------
## method = "length": alignment length / query length (used for the thesis
##                    genome-wide presence/absence matrix)
## method = "span":   (qend - qstart + 1) / query length
## Both gave identical pass/fail calls on the B. wexlerae DSM 19850 check.
add_query_coverage <- function(dt, method = c("length", "span")) {
  method <- match.arg(method)
  if (method == "length") dt[, qcov := 100 * length / qlen]
  else                    dt[, qcov := 100 * (qend - qstart + 1) / qlen]
  dt[]
}

## ---- Filtering and best hits -------------------------------------------
## Keep hits passing all three thresholds, then the single highest-bitscore
## hit for each genome-gene pair.
filter_best_hits <- function(dt, evalue_max = THRESHOLDS$evalue_max,
                             pident_min = THRESHOLDS$pident_min,
                             qcov_min = THRESHOLDS$qcov_min) {
  stopifnot(all(c("genome", "qseqid", "evalue", "pident", "qcov", "bitscore") %in% names(dt)))
  kept <- dt[evalue <= evalue_max & pident >= pident_min & qcov >= qcov_min]
  setorder(kept, genome, qseqid, -bitscore)
  kept[, .SD[1], by = .(genome, qseqid)]
}

## Binary genome x gene matrix; genomes/genes with no hits are kept as 0
presence_matrix <- function(best_hits, genomes, genes) {
  grid <- CJ(genome = unique(genomes), qseqid = unique(genes))
  grid[, present := as.integer(paste(genome, qseqid) %in%
                                 paste(best_hits$genome, best_hits$qseqid))]
  m <- dcast(grid, genome ~ qseqid, value.var = "present")
  mat <- as.matrix(m[, -1, with = FALSE])
  rownames(mat) <- m$genome
  mat
}

## Same shape, but cells hold % identity of the best hit (0 = absent)
identity_matrix <- function(best_hits, genomes, genes) {
  grid <- CJ(genome = unique(genomes), qseqid = unique(genes))
  grid <- merge(grid, best_hits[, .(genome, qseqid, pident)],
                by = c("genome", "qseqid"), all.x = TRUE)
  grid[is.na(pident), pident := 0]
  m <- dcast(grid, genome ~ qseqid, value.var = "pident")
  mat <- as.matrix(m[, -1, with = FALSE])
  rownames(mat) <- m$genome
  mat
}

## ---- Taxonomy (GTDB-Tk summary) ---------------------------------------
## Returns genome, genus, species from a file with columns user_genome and
## classification ("d__Bacteria;...;g__Blautia_A;s__Blautia_A wexlerae").
read_taxonomy <- function(path) {
  tax <- fread(path, sep = "\t", header = TRUE)
  stopifnot("Taxonomy file needs 'user_genome' and 'classification' columns" =
              all(c("user_genome", "classification") %in% names(tax)))
  out <- tax[, .(
    genome  = as.character(user_genome),
    genus   = str_remove(str_extract(classification, "g__[^;]*"), "^g__"),
    species = str_remove(str_extract(classification, "s__[^;]*"), "^s__")
  )]
  out[is.na(species) | species == "", species := "Unclassified"]
  out[is.na(genus) | genus == "", genus := "Unclassified"]
  unique(out, by = "genome")
}

## "Blautia_A wexlerae" -> "B. wexlerae"
short_species <- function(x) str_replace(x, "^Blautia(_[A-Za-z0-9]+)?\\s+", "B. ")

top_species <- function(taxonomy, n = 20) {
  taxonomy[species != "Unclassified", .N, by = species][order(-N)][seq_len(min(n, .N)), species]
}

## ---- Gene panel metadata -----------------------------------------------
read_gene_metadata <- function(path = "data/panel/gene_metadata.tsv") {
  fread(path, sep = "\t", header = TRUE, quote = "")
}

## ---- Species-level summaries ----------------------------------------------
## For each species x gene: number of genomes, number carrying the gene,
## fraction carrying it, and mean % identity (absent genomes count as 0).
species_gene_summary <- function(ident_mat, taxonomy) {
  long <- as.data.table(ident_mat, keep.rownames = "genome")
  long <- melt(long, id.vars = "genome", variable.name = "qseqid",
               value.name = "pident", variable.factor = FALSE)
  long <- merge(long, taxonomy[, .(genome, species)], by = "genome")
  long[, .(n_genomes = .N,
           n_hit = sum(pident > 0),
           presence_frac = mean(pident > 0),
           mean_pident = mean(pident)),
       by = .(species, qseqid)]
}

## Thesis tier definition, applied at species level: a gene is "present" in a
## species when >= 50% of its genomes carry it; across the top 20 species a
## gene is core in >= 18, accessory in 5-17 and rare in <= 4 species.
## For a different number of species the same proportions are used.
gene_tiers <- function(summary, presence_frac_min = 0.5, n_species = NULL) {
  if (is.null(n_species)) n_species <- uniqueN(summary$species)
  core_min <- ceiling(0.9 * n_species)
  acc_min  <- ceiling(0.25 * n_species)
  summary[, .(n_species_present = sum(presence_frac >= presence_frac_min)), by = qseqid][
    , tier := fifelse(n_species_present >= core_min, "core",
               fifelse(n_species_present >= acc_min, "accessory", "rare"))][]
}

## ---- Output helper -----------------------------------------------------------
ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)
