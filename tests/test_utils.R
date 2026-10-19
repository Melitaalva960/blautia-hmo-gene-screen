## Unit tests for R/utils.R.   Run from the repository root:
##   Rscript tests/test_utils.R
source("R/utils.R")

failures <- 0
check <- function(desc, ok) {
  cat(if (isTRUE(ok)) "PASS " else "FAIL ", desc, "\n", sep = "")
  if (!isTRUE(ok)) failures <<- failures + 1
}

tmp <- tempfile("tblastn"); dir.create(tmp)
## 13-column thesis format: qseqid qlen sseqid pident length mismatch gapopen
##                          qstart qend sstart send evalue bitscore
writeLines(c(
  "geneA\t300\tc1\t45\t300\t0\t0\t1\t300\t1\t900\t1e-50\t400",  # pass
  "geneA\t300\tc2\t40\t300\t0\t0\t1\t300\t1\t900\t1e-40\t300",  # pass, lower bitscore
  "geneB\t300\tc1\t60\t150\t0\t0\t1\t150\t1\t450\t1e-30\t200"   # fails coverage (50%)
), file.path(tmp, "G1_tblastn.tsv"))
writeLines(c(
  "geneA\t300\tc1\t25\t300\t0\t0\t1\t300\t1\t900\t1e-20\t150",  # fails identity
  "geneB\t300\tc9\t35\t280\t0\t0\t11\t290\t1\t840\t1e-3\t90",   # fails e-value
  "geneB\t300\tc1\t35\t280\t0\t0\t11\t290\t1\t840\t1e-10\t120"  # pass (93.3%)
), file.path(tmp, "G2_tblastn.tsv"))
invisible(file.create(file.path(tmp, "G3_tblastn.tsv")))                    # genome with no hits

hits <- read_tblastn_dir(tmp)
check("reads 13-column files and adds genome IDs", nrow(hits) == 6 && setequal(unique(hits$genome), c("G1", "G2")))
check("genome ID strips _tblastn suffix", genome_id_from_file("x/GUT_GENOME1_tblastn.tsv") == "GUT_GENOME1")

hits <- add_query_coverage(hits, "length")
check("coverage = length / qlen", abs(hits[genome == "G2" & qseqid == "geneB" & sseqid == "c1"]$qcov - 93.333) < 0.01)
span <- add_query_coverage(copy(hits), "span")
check("span coverage = (qend - qstart + 1) / qlen", abs(span[genome == "G2" & qseqid == "geneB" & sseqid == "c1"]$qcov - 93.333) < 0.01)

best <- filter_best_hits(hits)
check("keeps only passing genome-gene pairs",
      setequal(paste(best$genome, best$qseqid), c("G1 geneA", "G2 geneB")))
check("keeps the highest-bitscore hit", best[genome == "G1"]$bitscore == 400)

pa <- presence_matrix(best, c("G1", "G2", "G3"), c("geneA", "geneB", "geneC"))
check("matrix keeps genomes and genes without hits", all(dim(pa) == c(3, 3)))
check("matrix values are correct",
      pa["G1", "geneA"] == 1 && pa["G1", "geneB"] == 0 && sum(pa["G3", ]) == 0 && sum(pa[, "geneC"]) == 0)
idm <- identity_matrix(best, c("G1", "G2", "G3"), c("geneA", "geneB"))
check("identity matrix holds best-hit % identity", idm["G1", "geneA"] == 45 && idm["G3", "geneA"] == 0)

tax_file <- tempfile(fileext = ".tsv")
writeLines(c("user_genome\tclassification",
             "G1\td__Bacteria;g__Blautia_A;s__Blautia_A wexlerae",
             "G2\td__Bacteria;g__Blautia_A;s__Blautia_A wexlerae",
             "G3\td__Bacteria;g__Blautia;s__"), tax_file)
tax <- read_taxonomy(tax_file)
check("parses genus and species from GTDB classification",
      tax[genome == "G1"]$species == "Blautia_A wexlerae" && tax[genome == "G1"]$genus == "Blautia_A")
check("empty species becomes Unclassified", tax[genome == "G3"]$species == "Unclassified")
check("short species names", short_species("Blautia_A wexlerae") == "B. wexlerae" &&
        short_species("Blautia hansenii") == "B. hansenii")

summ <- species_gene_summary(idm, tax)
w <- summ[species == "Blautia_A wexlerae" & qseqid == "geneA"]
check("species summary: fraction and mean identity", w$presence_frac == 0.5 && w$mean_pident == 22.5)

## 20 species: gene in 19 -> core, in 10 -> accessory, in 2 -> rare
ts <- data.table(species = rep(paste0("s", 1:20), 3),
                 qseqid = rep(c("core_g", "acc_g", "rare_g"), each = 20),
                 presence_frac = c(rep(1, 19), 0, rep(1, 10), rep(0, 10), rep(1, 2), rep(0, 18)))
t <- gene_tiers(ts)
check("tiers follow the thesis definition (>=18 / 5-17 / <=4 of 20)",
      t[qseqid == "core_g"]$tier == "core" && t[qseqid == "acc_g"]$tier == "accessory" &&
        t[qseqid == "rare_g"]$tier == "rare")

meta <- read_gene_metadata()
check("gene panel metadata has 64 genes with unique labels", nrow(meta) == 64 && !anyDuplicated(meta$label))
groups <- fread("data/panel/gene_groups.tsv", sep = "\t", quote = "")
check("every grouped gene exists in the panel", all(groups$label %in% meta$label))

if (failures > 0) { cat(failures, "test(s) failed\n"); quit(status = 1) }
cat("All tests passed\n")
