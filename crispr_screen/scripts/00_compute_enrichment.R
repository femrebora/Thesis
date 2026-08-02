#!/usr/bin/env Rscript
# =============================================================================
# 00_compute_enrichment.R — compute the KEGG / Reactome / GO enrichment tables
# =============================================================================
# Produces the six result tables that figure scripts 10, 11 and 12 consume:
#
#   results/kegg_all_25.csv       results/kegg_Rpost_18.csv
#   results/reactome_all_25.csv   results/reactome_Rpost_18.csv
#   results/go_all_25.csv         results/go_Rpost_18.csv
#
# Scripts 10-12 deliberately read these precomputed tables so that figure
# rendering never depends on a live KEGG/Reactome API call. This script is the
# upstream step that creates them, and is the ONLY place enrichment statistics
# are computed.
#
# Inputs  : data/gene_lists/all_significant_25.txt
#           data/gene_lists/Rpost_vs_R0_significant_18.txt
# Requires: network access (KEGG REST API via clusterProfiler)
#
# Usage:
#   Rscript crispr_screen/scripts/00_compute_enrichment.R
#
# Note on the `Description_TR` column: the committed tables carry an extra
# Turkish-translation column used for the bilingual figure labels. It is added
# afterwards by build_enrichment_translations.py and is NOT produced here, so a
# fresh run of this script yields the same statistics with that column absent.
# Pass --keep-translations (default) to merge any existing Description_TR back in.
# =============================================================================

source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
})

GENE_SETS <- c(
  all_25   = "all_significant_25.txt",
  Rpost_18 = "Rpost_vs_R0_significant_18.txt"
)

# Enrichment parameters. pvalueCutoff = 0.1 with BH adjustment; no `universe`
# argument, so the background is each database's full annotated gene set.
# clusterProfiler returns every tested term in @result regardless of the cutoff,
# which is why the tables contain terms with p.adjust well above 0.1 — the
# figure scripts, not this script, decide what is displayed.
P_CUTOFF    <- 0.1
P_ADJ_METHOD <- "BH"

KEEP_TRANSLATIONS <- !("--no-keep-translations" %in% commandArgs(trailingOnly = TRUE))

# -- helpers -------------------------------------------------------------------
read_gene_list <- function(fname) {
  genes <- readLines(file.path(GENE_DIR, fname), warn = FALSE)
  genes <- trimws(genes)
  genes[nzchar(genes)]
}

to_entrez <- function(symbols) {
  map <- suppressWarnings(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  unmapped <- setdiff(symbols, map$SYMBOL)
  if (length(unmapped)) {
    warning("Unmapped symbols dropped: ", paste(unmapped, collapse = ", "),
            call. = FALSE)
  }
  map$ENTREZID
}

# Re-attach the Turkish description column from the previous table, matched on
# term ID, so regenerating the statistics never silently drops the translations.
merge_translations <- function(df, path) {
  if (!KEEP_TRANSLATIONS || !file.exists(path)) return(df)
  prev <- suppressWarnings(read_csv(path, show_col_types = FALSE))
  if (!"Description_TR" %in% names(prev)) return(df)
  df$Description_TR <- prev$Description_TR[match(df$ID, prev$ID)]
  df$Description_TR[is.na(df$Description_TR)] <- ""
  df
}

write_result <- function(obj, prefix, label) {
  if (is.null(obj) || nrow(obj@result) == 0) {
    message("  ", prefix, " / ", label, ": no terms returned — table not written")
    return(invisible(NULL))
  }
  res <- obj@result
  res$GeneSet <- label
  path <- file.path(TAB_DIR, paste0(prefix, "_", label, ".csv"))
  res <- merge_translations(res, path)
  dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
  write_csv(res, path)
  message("  ", basename(path), ": ", nrow(res), " terms")
  invisible(res)
}

# -- run -----------------------------------------------------------------------
for (label in names(GENE_SETS)) {
  symbols <- read_gene_list(GENE_SETS[[label]])
  entrez  <- to_entrez(symbols)
  message("Gene set '", label, "': ", length(symbols), " symbols -> ",
          length(entrez), " Entrez IDs")

  write_result(
    enrichKEGG(gene = entrez, organism = "hsa",
               pvalueCutoff = P_CUTOFF, pAdjustMethod = P_ADJ_METHOD),
    "kegg", label)

  write_result(
    enrichPathway(gene = entrez, organism = "human",
                  pvalueCutoff = P_CUTOFF, pAdjustMethod = P_ADJ_METHOD,
                  readable = TRUE),
    "reactome", label)

  write_result(
    enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "ALL",
             pvalueCutoff = P_CUTOFF, pAdjustMethod = P_ADJ_METHOD,
             readable = TRUE),
    "go", label)
}

message("Enrichment tables written to ", TAB_DIR)
