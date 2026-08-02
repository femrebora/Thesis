#!/usr/bin/env Rscript
# Module 5: Mutation & Genomic Alteration Landscape (TCGA-GBM)
# =================================================================
# Uses MC3 consensus MAF (pre-filtered to GBM cases) for mutation calls.
# CNA segment data mapped to gene-level via GDC if available.

# Auto-detect script directory for both Rscript and interactive modes
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
setwd(script_dir)

source("config.R")
source("utils.R")

library(maftools)
library(ggplot2)
library(reshape2)

dir.create(file.path(RESULTS_DIR, "05_mutation_landscape"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(RESULTS_DIR, "05_mutation_landscape")

# ---- 1. Load TCGA-GBM MAF data (pre-filtered from MC3 consensus) ----
message("[1/6] Loading TCGA-GBM mutation data (MC3 consensus calls)...")
maf_query <- read.delim(
  "gbm_mc3.maf",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

message("  Mutations loaded: ", nrow(maf_query), " variants in ",
        length(unique(maf_query$Tumor_Sample_Barcode)), " samples")

# ---- 2. Create MAF object ----
message("[2/6] Creating MAF object...")
maf <- read.maf(maf = maf_query, verbose = FALSE)

# Print summary
message("  MAF summary:")
message("    Samples: ", ncol(maf@data) - 1)  # approximate
message("    Variants: ", nrow(maf@data))

# ---- 3. Try CNA data download ----
message("[3/6] Attempting CNA data download...")
cna_available <- FALSE
cna_target <- NULL

tryCatch({
  library(TCGAbiolinks)

  # Try masked copy number segment data
  cna_query <- GDCquery(
    project = TCGA_PROJECT,
    data.category = "Copy Number Variation",
    data.type = "Masked Copy Number Segment",
    sample.type = "Primary Tumor"
  )

  if (nrow(cna_query$results) > 0) {
    GDCdownload(cna_query, files.per.chunk = 10)
    cna_data <- GDCprepare(cna_query)

    # Map segments to genes: compute mean segment value per gene
    # cna_data is a data.frame with columns: Chromosome, Start, End, Num_Probes, Segment_Mean, Sample
    message("  CNA segments loaded: ", nrow(cna_data))
    message("  CNA samples: ", length(unique(cna_data$Sample)))

    cna_available <- TRUE
  } else {
    message("  [NOTE] No CNA segment data available for GBM in GDC")
  }
}, error = function(e) {
  message("  [NOTE] CNA download failed: ", e$message)
  message("  [NOTE] Proceeding with mutation-only analysis")
})

# ---- 4. Determine alteration frequencies ----
message("[4/6] Computing alteration frequencies...")

# Get genes present in MAF
maf_genes_in_target <- intersect(ALL_GENES, unique(maf_query$Hugo_Symbol))
message("  Target genes with mutations: ",
        paste(maf_genes_in_target, collapse = ", "))

# Compile frequencies
alteration_summary <- data.frame(
  gene = ALL_GENES,
  class = GENE_CLASS[ALL_GENES],
  mutation_rate = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(alteration_summary))) {
  gene <- alteration_summary$gene[i]

  # Mutation frequency
  if (gene %in% unique(maf_query$Hugo_Symbol)) {
    gene_muts <- maf_query[maf_query$Hugo_Symbol == gene, ]
    n_mutated <- length(unique(gene_muts$Tumor_Sample_Barcode))
    n_total <- length(unique(maf_query$Tumor_Sample_Barcode))
    alteration_summary$mutation_rate[i] <- n_mutated / n_total
  } else {
    alteration_summary$mutation_rate[i] <- 0
  }
}

message("  Mutation frequencies:")
for (i in seq_len(nrow(alteration_summary))) {
  message(sprintf("    %s: %.1f%% (%s)",
    alteration_summary$gene[i],
    alteration_summary$mutation_rate[i] * 100,
    alteration_summary$class[i]))
}

# ---- 5. OncoPrint ----
message("[5/6] Generating OncoPrint...")

if (length(maf_genes_in_target) > 0) {
  # Interactive display
  oncoplot(
    maf = maf,
    genes = maf_genes_in_target,
    removeNonMutated = FALSE,
    titleText = "Mutation Landscape of CRISPR Screen Hits in TCGA-GBM",
    showTumorSampleBarcodes = FALSE,
    fontSize = 0.8
  )

  # Save PDF
  pdf(file.path(OUT, "oncoplot.pdf"), width = 10, height = 6)
  oncoplot(
    maf = maf,
    genes = maf_genes_in_target,
    removeNonMutated = FALSE,
    titleText = "Mutation Landscape of CRISPR Screen Hits in TCGA-GBM",
    showTumorSampleBarcodes = FALSE
  )
  dev.off()
  message("  Saved oncoplot PDF")

  # Save TIFF
  tiff(file.path(OUT, "oncoplot.tiff"), width = 10, height = 6,
       units = "in", res = 300, compression = "lzw")
  oncoplot(
    maf = maf,
    genes = maf_genes_in_target,
    removeNonMutated = FALSE,
    titleText = "Mutation Landscape of CRISPR Screen Hits in TCGA-GBM",
    showTumorSampleBarcodes = FALSE
  )
  dev.off()
  message("  Saved oncoplot TIFF")
} else {
  message("  [SKIP] No target genes found in MAF")
}

# ---- 6. Mutation frequency bar chart ----
message("[6/6] Generating mutation frequency summary figure...")

p_mut <- ggplot(alteration_summary,
  aes(x = reorder(gene, mutation_rate), y = mutation_rate * 100, fill = class)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Depleted" = "#2166AC", "Enriched" = "#B2182B",
                                "Other" = "gray60")) +
  labs(x = "", y = "Mutation Frequency (% of Samples)",
       fill = "CRISPR Class",
       title = "Somatic Mutation Frequency in TCGA-GBM",
       subtitle = paste0("MC3 consensus calls (",
                         length(unique(maf_query$Tumor_Sample_Barcode)),
                         " samples)")) +
  theme_gbm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(p_mut, file.path(OUT, "mutation_frequencies"), width = 10, height = 6)
message("  Saved mutation frequency figure")

# ---- 7. Co-occurrence with GBM drivers ----
message("  Computing co-occurrence with GBM drivers...")
drivers <- c("TP53", "PTEN", "EGFR", "IDH1", "NF1", "RB1")
drivers_in_maf <- intersect(drivers, unique(maf_query$Hugo_Symbol))

if (length(maf_genes_in_target) > 0 && length(drivers_in_maf) > 0) {
  message("  Driver genes found: ", paste(drivers_in_maf, collapse = ", "))

  cooccur <- somaticInteractions(
    maf = maf,
    top = min(30, length(unique(maf_query$Hugo_Symbol))),
    genes = c(maf_genes_in_target, drivers_in_maf)
  )

  pdf(file.path(OUT, "cooccurrence.pdf"), width = 10, height = 8)
  plot(cooccur)
  dev.off()
  message("  Saved co-occurrence PDF")

  tiff(file.path(OUT, "cooccurrence.tiff"), width = 10, height = 8,
       units = "in", res = 300, compression = "lzw")
  plot(cooccur)
  dev.off()
  message("  Saved co-occurrence TIFF")
} else {
  message("  [SKIP] No target genes or drivers for co-occurrence analysis")
}

# ---- 8. Lollipop plots for top mutated genes ----
message("  Generating lollipop plots for top mutated target genes...")
if (length(maf_genes_in_target) > 0) {
  # Sort by mutation count
  gene_mut_counts <- sapply(maf_genes_in_target, function(g) {
    sum(maf_query$Hugo_Symbol == g)
  })
  gene_mut_counts <- sort(gene_mut_counts, decreasing = TRUE)

  for (gene_name in names(gene_mut_counts)[1:min(3, length(gene_mut_counts))]) {
    tryCatch({
      pdf(file.path(OUT, paste0("lollipop_", gene_name, ".pdf")), width = 10, height = 4)
      lollipopPlot(
        maf = maf,
        gene = gene_name,
        showMutationRate = TRUE,
        labelPos = "all"
      )
      dev.off()
      message(sprintf("    Lollipop: %s (%d mutations)", gene_name, gene_mut_counts[gene_name]))
    }, error = function(e) {
      message(sprintf("    [SKIP] Lollipop for %s: %s", gene_name, e$message))
    })
  }
}

# ---- 9: Save tables ----
save_table(alteration_summary, file.path(OUT, "alteration_table"))
message("  Saved alteration frequency table")

# ---- Summary ----
message("\n==== Module 5 Summary ====")
message("Total samples: ", length(unique(maf_query$Tumor_Sample_Barcode)))
message("Genes with mutations: ", length(maf_genes_in_target), "/", length(ALL_GENES))
message("Genes mutated: ", paste(maf_genes_in_target, collapse = ", "))

message("[done] Module 5 complete. Outputs in ", OUT)
