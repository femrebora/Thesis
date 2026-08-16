#!/usr/bin/env Rscript
# Entry point for the CRISPR/Cas9 metabolic screen figure pipeline.
# Runs each figure script in order, records pass/fail and timing, and writes a
# sessionInfo() snapshot for reproducibility.
#
# Usage (from anywhere inside the repo):
#   Rscript crispr_screen/scripts/run_crispr_analysis.R
#
# Inputs (data/, git-ignored — see crispr_screen/data/README.md):
#   - count_table.txt                 MAGeCK count matrix (sgRNA x sample)
#   - <comparison>.gene_summary.tsv   MAGeCK gene-level test output (x3)
#   - gene_lists/*.txt                gene symbols fed to enrichment (committed)

suppressPackageStartupMessages(library(here))

module    <- here::here("crispr_screen")
scriptdir <- file.path(module, "scripts")
datadir   <- file.path(module, "data")

# Fail early and clearly if the (restricted) input data is not present.
required <- c(
  file.path(datadir, "count_table.txt"),
  file.path(datadir, "Rpost_vs_R0.gene_summary.tsv"),
  file.path(datadir, "Spost_vs_Rpost.gene_summary.tsv"),
  file.path(datadir, "Spost_vs_S0.gene_summary.tsv")
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing required input(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nThese files are restricted and are not distributed with the public ",
       "repository. For schemas and placement, see crispr_screen/data/README.md. ",
       "Module overview: crispr_screen/README.md. ",
       "A fresh public clone can inspect code, committed gene lists, enrichment ",
       "tables, and baseline figures, but cannot regenerate CRISPR thesis figures ",
       "without these inputs.", call. = FALSE)
}

# Final-thesis figure scripts, in dependency order. 05b (palette test) and
# 07_replicate_corr are intentionally excluded — neither is a final thesis figure.
scripts <- c(
  "01_volcano.R", "02_rank_lfc.R", "03_heatmap.R", "04_barplot.R",
  "05_qc_metrics.R", "06_sgrna_profiles.R", "08_cross_comparison.R",
  "09_summary_table.R", "10_kegg_enrichment.R", "11_reactome_enrichment.R",
  "12_go_enrichment.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
results <- data.frame(script = scripts, status = NA_character_, seconds = NA_real_)

for (i in seq_along(scripts)) {
  s <- scripts[i]
  message("=== Running ", s, " ===")
  t0 <- Sys.time()
  code <- system2(rscript, file.path(scriptdir, s), stdout = "", stderr = "")
  results$seconds[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  results$status[i]  <- if (code == 0) "PASS" else paste0("FAIL(", code, ")")
}

cat("\n================ SUMMARY ================\n")
print(results, row.names = FALSE)
cat("=========================================\n")

si_path <- file.path(module, "results",
                     paste0("sessionInfo", Sys.getenv("THESIS_SIG_SUFFIX", ""), ".txt"))
dir.create(dirname(si_path), showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), si_path)
message("Wrote ", si_path)

if (any(grepl("FAIL", results$status))) quit(status = 1)
