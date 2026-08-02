#!/usr/bin/env Rscript
# Module 4: Immune Infiltration Correlation Analysis
# =================================================================
# Implements a manual MCP-counter-like approach using well-established
# immune cell-type marker gene signatures. Since immunedeconv could not
# be installed (GitHub API unreachable), we compute per-sample enrichment
# scores as the mean expression of each cell type's marker genes,
# z-scored across samples.
#
# Expression data is loaded directly from STAR-counts TSV files
# (previously downloaded by Module 1 via GDCdownload), avoiding any
# network dependency.
#
# Marker gene sets from: Becht et al. 2016 (MCP-counter) and
# Bindea et al. 2013 (CIBERSORT-like signatures).
#
# Outputs:
#   immune_infiltration_scores.csv       - Per-sample immune cell scores
#   immune_correlation_heatmap.{pdf,tiff} - Gene-immune correlation heatmap
#   immune_scatter_top.{pdf,tiff}        - Top 9 significant correlations
#   immune_corr_table.csv                - Full correlation table with FDR

# Set working directory to script location
script_dir <- tryCatch(
  dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])),
  error = function(e) getwd()
)
if (is.na(script_dir) || script_dir == "") script_dir <- getwd()
setwd(script_dir)

source("config.R")
source("utils.R")

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(data.table)

dir.create(file.path(RESULTS_DIR, "04_immune_infiltration"), showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(RESULTS_DIR, "04_immune_infiltration")

# ---- Marker gene definitions for immune cell types ----
# These are the well-established MCP-counter markers (Becht et al. 2016, Genome Biology)
# supplemented with additional validated markers from Bindea et al. 2013.
IMMUNE_MARKERS <- list(
  `T cells`               = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD2", "CD5", "CD7"),
  `CD8 T cells`           = c("CD8A", "CD8B", "GZMK", "GZMA", "PRF1"),
  `NK cells`              = c("NKG7", "KLRD1", "KLRB1", "NCR1", "GNLY", "KIR2DL1", "KIR2DL3"),
  `B cells`               = c("CD19", "CD79A", "CD79B", "MS4A1", "IGHM", "IGHD"),
  `Monocytes`             = c("CD14", "FCGR3A", "CSF1R", "LYZ", "CTSS"),
  `Macrophages`           = c("CD68", "CD163", "MSR1", "C1QA", "C1QB", "MARCO"),
  `Neutrophils`           = c("FCGR3B", "CEACAM3", "S100A8", "S100A9", "CSF3R"),
  `Dendritic cells`       = c("FCER1A", "CLEC10A", "CLEC9A", "CD1C", "LAMP3", "BATF3"),
  `Mast cells`            = c("KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  `Eosinophils`           = c("PRG2", "PRG3", "CLC", "CCR3", "IL5RA"),
  `Fibroblasts`           = c("COL1A1", "COL1A2", "FAP", "PDGFRA", "PDGFRB", "DCN", "LUM"),
  `Endothelial cells`     = c("PECAM1", "VWF", "CDH5", "ENG", "KDR", "FLT1")
)

# ---- 1. Load TCGA-GBM expression from TSV files ----
message("[1/5] Loading TCGA-GBM expression from STAR-counts TSV files...")

# Find all TSV files from the GDC download
tsv_dir <- file.path("GDCdata", TCGA_PROJECT, "Transcriptome_Profiling",
                     "Gene_Expression_Quantification")
tsv_files <- list.files(tsv_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
message("  Found ", length(tsv_files), " TSV files")

if (length(tsv_files) == 0) {
  stop("No TSV files found. Run Module 1 first to download TCGA-GBM data.")
}

# Read the first file to get gene annotations
first_file <- tsv_files[1]
header <- readLines(first_file, n = 1)
dt_first <- fread(first_file, skip = 1, sep = "\t", header = TRUE)
gene_ids   <- dt_first$gene_id
gene_names <- dt_first$gene_name

# Extract sample IDs from file paths
# GDC file names contain UUIDs; we extract the TCGA barcode from the file manifest
# The sample ID is embedded in the filename pattern or the parent directory UUID
# We'll use a UUID-to-sample mapping approach

# Build expression matrix from all files (using tpm_unstranded, column 7)
message("  Building expression matrix (tpm_unstranded)...")

# Strategy: read all files, extract tpm_unstranded column
# First, collect all gene IDs for row names
all_genes <- dt_first$gene_id
n_genes <- length(all_genes)
n_samples <- length(tsv_files)

tpm_matrix <- matrix(NA_real_, nrow = n_genes, ncol = n_samples)
rownames(tpm_matrix) <- all_genes

# We need sample IDs. GDC file naming doesn't embed TCGA barcode directly.
# The parent directories are UUIDs. We need to get TCGA barcodes from somewhere.
# Option: the file path contains the case UUID; we can use that.
# But for a simpler approach, we'll just use sample indices.
# Better approach: read all files, use the manifest UUIDs as colnames,
# and then when matching against our gene list, we only need TPM values.

col_names <- basename(dirname(tsv_files))  # Use UUID directory names
message("  Loading ", n_samples, " samples x ", n_genes, " genes...")

pb_tick <- ceiling(n_samples / 20)
for (i in seq_along(tsv_files)) {
  dt <- fread(tsv_files[i], skip = 1, sep = "\t", header = TRUE,
              select = c("gene_id", "tpm_unstranded"))
  # Match and fill
  idx <- match(all_genes, dt$gene_id)
  tpm_matrix[, i] <- dt$tpm_unstrand[idx]
  if (i %% pb_tick == 0) message("    ", i, "/", n_samples, " files processed")
}
colnames(tpm_matrix) <- col_names

message("  Expression matrix: ", nrow(tpm_matrix), " genes x ", ncol(tpm_matrix), " samples")

# ---- 2. Map to gene symbols and deduplicate ----
message("[2/5] Mapping to gene symbols...")

# Filter to protein-coding and known genes (exclude N_*, __no_feature, etc.)
keep <- !grepl("^_", all_genes) & !grepl("^N_", all_genes)
tpm_matrix <- tpm_matrix[keep, ]
gene_names <- gene_names[keep]

# Map to symbol, average duplicates
gene_symbols <- as.character(gene_names)

duplicated_symbols <- unique(gene_symbols[duplicated(gene_symbols)])
if (length(duplicated_symbols) > 0) {
  message("  Averaging ", length(duplicated_symbols), " duplicate gene symbols...")
}

# Sum by symbol, then divide by group size to get mean
tpm_dedup <- rowsum(tpm_matrix, gene_symbols, reorder = TRUE)
group_sizes <- as.matrix(rowsum(rep(1, length(gene_symbols)), gene_symbols, reorder = TRUE))
tpm_dedup <- tpm_dedup / as.numeric(group_sizes[, 1])

message("  Deduplicated matrix: ", nrow(tpm_dedup), " genes x ", ncol(tpm_dedup), " samples")

# ---- 3. Compute immune infiltration scores (manual MCP-counter) ----
message("[3/5] Computing immune cell enrichment scores...")

# Log2 transform TPM
mat_log2 <- log2(tpm_dedup + 1)

# For each immune cell type, compute mean expression of its marker genes
# and then z-score across samples (as MCP-counter does)
immune_scores <- matrix(NA_real_, nrow = length(IMMUNE_MARKERS), ncol = ncol(mat_log2))
rownames(immune_scores) <- names(IMMUNE_MARKERS)
colnames(immune_scores) <- colnames(mat_log2)

all_found <- character()
all_missing <- character()

for (ct in names(IMMUNE_MARKERS)) {
  markers <- IMMUNE_MARKERS[[ct]]
  found <- intersect(markers, rownames(mat_log2))
  missing <- setdiff(markers, rownames(mat_log2))

  all_found <- c(all_found, found)
  all_missing <- c(all_missing, missing)

  if (length(found) >= 3) {
    # Mean expression across marker genes
    ct_score <- colMeans(mat_log2[found, , drop = FALSE], na.rm = TRUE)
  } else if (length(found) > 0) {
    message("  NOTE: ", ct, " only has ", length(found), "/", length(markers), " markers")
    ct_score <- colMeans(mat_log2[found, , drop = FALSE], na.rm = TRUE)
  } else {
    message("  WARNING: ", ct, " has NO markers in expression data")
    next
  }

  # Z-score normalization across samples
  if (sd(ct_score, na.rm = TRUE) > 0) {
    immune_scores[ct, ] <- (ct_score - mean(ct_score, na.rm = TRUE)) / sd(ct_score, na.rm = TRUE)
  } else {
    # If all identical (rare), skip this cell type
    message("  NOTE: ", ct, " scores have zero variance, skipping")
  }
}

# Remove cell types that ended up with no valid scores
valid_ct <- rowSums(!is.na(immune_scores)) > 0
if (sum(valid_ct) < nrow(immune_scores)) {
  removed_ct <- rownames(immune_scores)[!valid_ct]
  message("  Removing cell types with invalid scores: ", paste(removed_ct, collapse = ", "))
  immune_scores <- immune_scores[valid_ct, , drop = FALSE]
}

message("  Marker genes found: ", length(unique(all_found)), " / ",
        length(unique(c(all_found, all_missing))))
message("  Immune scores matrix: ", nrow(immune_scores), " cell types x ", ncol(immune_scores), " samples")

# Save immune scores to CSV
scores_for_csv <- as.data.frame(t(immune_scores))
scores_for_csv <- cbind(sample_id = rownames(scores_for_csv), scores_for_csv)
write.csv(scores_for_csv, file.path(OUT, "immune_infiltration_scores.csv"), row.names = FALSE)
message("  Saved: immune_infiltration_scores.csv")

# ---- 4. Compute correlation with gene expression ----
message("[4/5] Computing gene-immune correlations...")

# Extract expression for our target genes
common_genes <- intersect(ALL_GENES, rownames(tpm_dedup))
message("  Target genes available: ", paste(common_genes, collapse = ", "))

expr_target <- mat_log2[common_genes, , drop = FALSE]

immune_cell_types <- rownames(immune_scores)

cor_matrix <- matrix(NA_real_, nrow = length(common_genes), ncol = length(immune_cell_types))
rownames(cor_matrix) <- common_genes
colnames(cor_matrix) <- immune_cell_types
pval_matrix <- cor_matrix

for (g in common_genes) {
  for (ct in immune_cell_types) {
    common_samples <- intersect(colnames(expr_target), colnames(immune_scores))
    if (length(common_samples) >= 5) {
      g_vec <- as.numeric(expr_target[g, common_samples])
      ct_vec <- as.numeric(immune_scores[ct, common_samples])
      valid <- !is.na(g_vec) & !is.na(ct_vec)
      if (sum(valid) >= 5) {
        test <- cor.test(g_vec[valid], ct_vec[valid], method = "spearman")
        cor_matrix[g, ct] <- test$estimate
        pval_matrix[g, ct] <- test$p.value
      }
    }
  }
}

# FDR correction across all tests
pval_adj <- matrix(
  p.adjust(as.vector(pval_matrix), method = "BH"),
  nrow = nrow(pval_matrix),
  dimnames = dimnames(pval_matrix)
)

message("  Computed ", sum(!is.na(cor_matrix)), " correlations")
message("  Significant (FDR < 0.05): ", sum(pval_adj < 0.05, na.rm = TRUE))

# ---- 5. Generate outputs ----
message("[5/5] Generating figures and tables...")

# ----- 5a. Correlation heatmap -----
# Remove any rows/columns that are all NA
valid_rows <- rowSums(is.na(cor_matrix)) < ncol(cor_matrix)
valid_cols <- colSums(is.na(cor_matrix)) < nrow(cor_matrix)
cor_plot <- cor_matrix[valid_rows, valid_cols, drop = FALSE]
pval_plot <- pval_matrix[valid_rows, valid_cols, drop = FALSE]
pval_adj_plot <- pval_adj[valid_rows, valid_cols, drop = FALSE]

if (nrow(cor_plot) > 0 && ncol(cor_plot) > 0) {
  # Row annotation for gene class
  ha_row <- rowAnnotation(
    Class = GENE_CLASS[rownames(cor_plot)],
    col = list(Class = c("Depleted" = "#377EB8", "Enriched" = "#E41A1C", "Other" = "#999999")),
    show_annotation_name = TRUE
  )

  # Define colour range with symmetric limits
  max_abs <- max(abs(cor_plot), na.rm = TRUE)
  max_abs <- max(max_abs, 0.1)  # ensure minimum range
  col_breaks <- c(-max_abs, 0, max_abs)

  hmap <- Heatmap(
    cor_plot,
    name = "Spearman rho",
    col = colorRamp2(col_breaks, c("#2166AC", "white", "#B2182B")),
    right_annotation = ha_row,
    row_names_side = "left",
    column_names_rot = 45,
    column_title = "Gene-Immune Cell Infiltration Correlation (TCGA-GBM)",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(pval_adj_plot[i, j]) && pval_adj_plot[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 14, fontface = "bold"))
      }
    },
    heatmap_legend_param = list(title = "Spearman rho", direction = "horizontal")
  )

  pdf(file.path(OUT, "immune_correlation_heatmap.pdf"), width = 12, height = 8)
  draw(hmap, heatmap_legend_side = "bottom")
  dev.off()

  tiff(file.path(OUT, "immune_correlation_heatmap.tiff"), width = 12, height = 8,
       units = "in", res = 300, compression = "lzw")
  draw(hmap, heatmap_legend_side = "bottom")
  dev.off()

  message("  Saved: immune_correlation_heatmap.{pdf,tiff}")
} else {
  message("  WARNING: No valid correlations to plot in heatmap")
}

# ----- 5b. Top correlation scatter plots -----
cor_long <- data.frame(
  gene = rep(rownames(cor_matrix), ncol(cor_matrix)),
  cell_type = rep(colnames(cor_matrix), each = nrow(cor_matrix)),
  rho = as.vector(cor_matrix),
  p_adj = as.vector(pval_adj),
  stringsAsFactors = FALSE
)
cor_long <- cor_long[!is.na(cor_long$rho) & !is.na(cor_long$p_adj), ]
cor_long <- cor_long[order(-abs(cor_long$rho)), ]
top_hits <- head(cor_long[cor_long$p_adj < 0.05, ], 9)

if (nrow(top_hits) > 0) {
  message("  Generating ", nrow(top_hits), " scatter plots...")
  scatter_plots <- list()
  for (i in 1:nrow(top_hits)) {
    g   <- top_hits$gene[i]
    ct  <- top_hits$cell_type[i]
    common_samples <- intersect(colnames(expr_target), colnames(immune_scores))

    df <- data.frame(
      expression   = as.numeric(expr_target[g, common_samples]),
      infiltration = as.numeric(immune_scores[ct, common_samples])
    )
    df <- df[complete.cases(df), ]

    scatter_plots[[i]] <- ggplot(df, aes(x = expression, y = infiltration)) +
      geom_point(alpha = 0.4, size = 1.5, color = "#377EB8") +
      geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", linewidth = 1) +
      labs(
        x = paste0(g, " log2(TPM+1)"),
        y = paste0(ct, " score (z)"),
        title = sprintf("%s vs %s: rho=%.2f, FDR=%.3f", g, ct, top_hits$rho[i], top_hits$p_adj[i])
      ) +
      theme_gbm(base_size = 9)
  }

  n_col <- min(3, length(scatter_plots))
  combo <- cowplot::plot_grid(plotlist = scatter_plots, ncol = n_col)
  save_figure(combo, file.path(OUT, "immune_scatter_top"),
              width = 5 * n_col,
              height = 4 * ceiling(length(scatter_plots) / n_col))
  message("  Saved: immune_scatter_top.{pdf,tiff}")
} else {
  message("  No significant correlations (FDR < 0.05). Using top 9 by absolute rho instead.")
  top_hits <- head(cor_long, 9)
  if (nrow(top_hits) > 0) {
    scatter_plots <- list()
    for (i in 1:nrow(top_hits)) {
      g   <- top_hits$gene[i]
      ct  <- top_hits$cell_type[i]
      common_samples <- intersect(colnames(expr_target), colnames(immune_scores))

      df <- data.frame(
        expression   = as.numeric(expr_target[g, common_samples]),
        infiltration = as.numeric(immune_scores[ct, common_samples])
      )
      df <- df[complete.cases(df), ]

      scatter_plots[[i]] <- ggplot(df, aes(x = expression, y = infiltration)) +
        geom_point(alpha = 0.4, size = 1.5, color = "#377EB8") +
        geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", linewidth = 1) +
        labs(
          x = paste0(g, " log2(TPM+1)"),
          y = paste0(ct, " score (z)"),
          title = sprintf("%s vs %s: rho=%.2f, p=%.3f", g, ct, top_hits$rho[i], top_hits$p_adj[i])
        ) +
        theme_gbm(base_size = 9)
    }

    n_col <- min(3, length(scatter_plots))
    combo <- cowplot::plot_grid(plotlist = scatter_plots, ncol = n_col)
    save_figure(combo, file.path(OUT, "immune_scatter_top"),
                width = 5 * n_col,
                height = 4 * ceiling(length(scatter_plots) / n_col))
    message("  Saved: immune_scatter_top.{pdf,tiff} (using top correlations)")
  } else {
    message("  No correlations available at all. Skipping scatter plots.")
  }
}

# ----- 5c. Save correlation table -----
cor_table <- data.frame(
  gene = rep(rownames(cor_matrix), ncol(cor_matrix)),
  cell_type = rep(colnames(cor_matrix), each = nrow(cor_matrix)),
  spearman_rho = round(as.vector(cor_matrix), 3),
  p_value = format_pval(as.vector(pval_matrix)),
  fdr_adj = format_pval(as.vector(pval_adj)),
  gene_class = GENE_CLASS[rep(rownames(cor_matrix), ncol(cor_matrix))],
  stringsAsFactors = FALSE
)
cor_table <- cor_table[order(abs(cor_table$spearman_rho), decreasing = TRUE), ]

save_table(cor_table, file.path(OUT, "immune_corr_table"))
message("  Saved: immune_corr_table.csv")

message("[done] Module 4 complete. Outputs in ", OUT)
