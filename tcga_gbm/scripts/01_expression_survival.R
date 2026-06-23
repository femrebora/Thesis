#!/usr/bin/env Rscript
# ============================================================================
# Module 01: TCGA-GBM Expression & Kaplan-Meier Survival Analysis
# ============================================================================
#
# Purpose:
#   Download and process TCGA-GBM RNA-seq expression data via TCGAbiolinks,
#   generate publication-ready expression visualizations, run Kaplan-Meier
#   survival analysis for target metabolic genes from CRISPR screen, and
#   produce a Cox regression forest plot.
#
# Target genes: 12 metabolic genes (4 depleted, 5 enriched, 3 other)
#   from TRAIL-resistance CRISPR/Cas9 screen (Cingoz et al., 2021).
#
# Analysis:   Tumor-only expression (TCGA-GBM has insufficient normal samples)
# Survival:   Overall survival, median-split high vs low expression
#
# Output:     figures/ (PDF + PNG 600dpi + TIFF 600dpi)
#             tables/  (CSV + TSV + XLSX)
#             01_expression_survival_summary.md
#             sessionInfo.txt
#
# Usage:      Rscript 01_expression_survival.R
#             (must be run from analysis/ directory)
#
# Version:    2.0
# ============================================================================

# ---- 0. Setup ----
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))),
  error = function(e) getwd()
)
if (is.na(SCRIPT_DIR) || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
setwd(SCRIPT_DIR)

source("config.R")
source("utils.R")

# ============================================================================
# 0. Package loading
# ============================================================================

load_packages <- function() {
  cli::cli_h1("Module 01: TCGA-GBM Expression & Survival Analysis")
  cli::cli_text("Checking required packages...")

  cran_pkgs <- c(
    "ggplot2", "dplyr", "tidyr", "tibble", "readr", "stringr", "purrr",
    "survival", "survminer", "maxstat",
    "ggpubr", "rstatix", "ggbeeswarm", "ggrepel", "patchwork",
    "broom", "viridis", "scales",
    "openxlsx", "cli", "glue", "here"
  )
  bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "ComplexHeatmap")

  # CRAN
  for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_alert_info("Installing {.pkg {pkg}}...")
      install.packages(pkg, repos = "https://cran.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }

  # Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.r-project.org")
  }
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_alert_info("Installing Bioconductor package {.pkg {pkg}}...")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    library(pkg, character.only = TRUE)
  }

  cli::cli_alert_success("All packages loaded")
  invisible(TRUE)
}

# ============================================================================
# 1. Data acquisition
# ============================================================================

download_tcga_expression <- function(project = TCGA_PROJECT, cache_dir = CACHE_DIR) {
  cli::cli_h2("Downloading TCGA-GBM RNA-seq data")

  query_exp <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  cli::cli_alert_info("Query returned {.val {nrow(query_exp$results[[1]])}} files")

  # Download with retry
  download_success <- FALSE
  for (chunk_size in c(10, 5, 3)) {
    tryCatch({
      GDCdownload(query_exp, files.per.chunk = chunk_size)
      download_success <- TRUE
      cli::cli_alert_success("Downloaded with files.per.chunk = {.val {chunk_size}}")
      break
    }, error = function(e) {
      cli::cli_alert_warning("Download with chunk={.val {chunk_size}} failed: {e$message}")
    })
  }
  if (!download_success) stop("Failed to download TCGA expression data after retries")

  exp_data <- GDCprepare(query_exp, summarizedExperiment = TRUE)
  cli::cli_alert_success("GDCprepare complete: {.val {ncol(exp_data)}} samples, {.val {nrow(exp_data)}} features")
  exp_data
}

download_tcga_clinical <- function(project = TCGA_PROJECT) {
  cli::cli_h2("Downloading TCGA-GBM clinical data")

  clinical <- GDCquery_clinic(project = project, type = "clinical")
  cli::cli_alert_info("Clinical data: {.val {nrow(clinical)}} patients, {.val {ncol(clinical)}} fields")

  required_cols <- c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up")
  missing_cols <- setdiff(required_cols, colnames(clinical))
  if (length(missing_cols) > 0) {
    cli::cli_alert_warning("Missing clinical columns: {.val {missing_cols}}")
  }

  clinical
}

# ============================================================================
# 2. Expression matrix preparation
# ============================================================================

prepare_expression_matrix <- function(se) {
  cli::cli_h2("Preparing expression matrix")

  # Detect available assay
  available <- names(assays(se))
  cli::cli_alert_info("Available assays: {.val {available}}")

  # Try TPM assays in order of preference
  preferred_assays <- c("tpm_unstranded", "tpm_unstrand")
  assay_used <- NULL

  for (a in preferred_assays) {
    if (a %in% available) {
      assay_used <- a
      break
    }
  }

  if (is.null(assay_used)) {
    if ("unstranded" %in% available) {
      cli::cli_alert_warning("TPM assay not found. Using raw 'unstranded' counts — will convert to log2(CPM+1)")
      assay_used <- "unstranded"
    } else {
      # Use first available as last resort
      assay_used <- available[1]
      cli::cli_alert_warning("No standard assay found. Using {.val {assay_used}}")
    }
  }
  cli::cli_alert_success("Assay selected: {.val {assay_used}}")

  expr_raw <- assay(se, assay_used)
  gene_annot <- rowData(se)

  # Map gene symbols
  raw_genes <- gene_annot$gene_name[match(rownames(expr_raw), rownames(gene_annot))]
  rownames(expr_raw) <- raw_genes

  # Count NA genes
  na_genes <- sum(is.na(raw_genes) | raw_genes == "")
  if (na_genes > 0) {
    cli::cli_alert_info("Removing {.val {na_genes}} features with missing gene symbols")
    expr_raw <- expr_raw[!is.na(raw_genes) & raw_genes != "", , drop = FALSE]
    raw_genes <- raw_genes[!is.na(raw_genes) & raw_genes != ""]
    rownames(expr_raw) <- raw_genes
  }

  # Handle duplicate gene symbols: keep row with highest mean expression
  dup_genes <- names(which(table(rownames(expr_raw)) > 1))
  if (length(dup_genes) > 0) {
    cli::cli_alert_info("Resolving {.val {length(dup_genes)}} duplicated gene symbol{?s}...")
    genes_kept <- c()
    genes_dropped <- c()

    for (g in dup_genes) {
      rows <- which(rownames(expr_raw) == g)
      row_means <- rowMeans(expr_raw[rows, , drop = FALSE], na.rm = TRUE)
      keep_idx <- rows[which.max(row_means)]
      genes_kept <- c(genes_kept, keep_idx)
      genes_dropped <- c(genes_dropped, setdiff(rows, keep_idx))
    }

    expr_raw <- expr_raw[-genes_dropped, , drop = FALSE]
    cli::cli_alert_info("  Kept {.val {length(dup_genes)}} highest-mean rows, dropped {.val {length(genes_dropped)}} duplicates")
  }

  # Log2 transform: TPM + 1
  expr_log2 <- log2(expr_raw + 1)
  cli::cli_alert_success("Expression matrix: {.val {nrow(expr_log2)}} genes × {.val {ncol(expr_log2)}} samples")
  cli::cli_alert_info("Values are log2({assay_used} + 1)")

  expr_log2
}

detect_sample_types <- function(se) {
  cli::cli_h2("Detecting TCGA sample types")

  # Extract TCGA barcodes
  barcodes <- NULL
  if ("barcode" %in% colnames(colData(se))) {
    barcodes <- colData(se)$barcode
  } else if ("sample" %in% colnames(colData(se))) {
    barcodes <- colData(se)$sample
  } else if ("aliquot_submitter_id" %in% colnames(colData(se))) {
    barcodes <- colData(se)$aliquot_submitter_id
  } else {
    barcodes <- colnames(se)
  }

  sample_info <- parse_tcga_sample_type(barcodes)
  counts <- table(sample_info$type)

  cli::cli_alert_info("Sample types detected:")
  for (t in names(counts)) {
    cli::cli_li("{t}: {.val {counts[t]}}")
  }

  # Warn about normal sample insufficiency
  n_normal <- sum(counts[names(counts) == "Solid Tissue Normal"])
  if (n_normal < MIN_NORMAL_SAMPLES) {
    cli::cli_alert_warning(
      "Only {.val {n_normal}} normal sample{?s} found (threshold: {.val {MIN_NORMAL_SAMPLES}}). ",
      "Tumor-normal comparison will be skipped."
    )
  }

  list(
    barcodes = barcodes,
    sample_code = sample_info$code,
    sample_type = sample_info$type,
    counts = as.data.frame(counts)
  )
}

select_and_average_tumor_samples <- function(expr_matrix, sample_info) {
  cli::cli_h2("Selecting primary tumor samples and averaging by patient")

  # Keep only Primary Tumor
  tumor_idx <- which(sample_info$sample_type == "Primary Tumor")
  cli::cli_alert_info("Primary Tumor samples: {.val {length(tumor_idx)}} / {.val {ncol(expr_matrix)}}")

  if (length(tumor_idx) == 0) {
    stop("No Primary Tumor samples found. Check TCGA barcode parsing.")
  }

  expr_tumor <- expr_matrix[, tumor_idx, drop = FALSE]
  tumor_barcodes <- sample_info$barcodes[tumor_idx]

  # Average by patient
  patient_ids <- parse_tcga_patient(tumor_barcodes)
  unique_patients <- unique(patient_ids)

  cli::cli_alert_info("Averaging {.val {length(patient_ids)}} aliquots → {.val {length(unique_patients)}} unique patients")

  expr_patient <- matrix(
    NA_real_,
    nrow = nrow(expr_tumor),
    ncol = length(unique_patients),
    dimnames = list(rownames(expr_tumor), unique_patients)
  )

  for (i in seq_along(unique_patients)) {
    pid <- unique_patients[i]
    cols <- which(patient_ids == pid)
    if (length(cols) == 1) {
      expr_patient[, i] <- expr_tumor[, cols]
    } else {
      expr_patient[, i] <- rowMeans(expr_tumor[, cols, drop = FALSE], na.rm = TRUE)
    }
  }

  cli::cli_alert_success("Patient-level expression: {.val {ncol(expr_patient)}} patients × {.val {nrow(expr_patient)}} genes")
  expr_patient
}

extract_target_genes <- function(expr_matrix, target_genes = ALL_GENES) {
  cli::cli_h2("Extracting target genes")

  found_genes <- intersect(target_genes, rownames(expr_matrix))
  missing_genes <- setdiff(target_genes, rownames(expr_matrix))

  cli::cli_alert_info("Target genes: {.val {length(target_genes)}}")
  cli::cli_alert_success("Found: {.val {length(found_genes)}} — {.val {found_genes}}")

  if (length(missing_genes) > 0) {
    cli::cli_alert_warning("Missing: {.val {length(missing_genes)}} — {.val {missing_genes}}")
    cli::cli_li("Reason: gene symbol not present in TCGA STAR-Counts annotation")
  }

  expr_target <- expr_matrix[found_genes, , drop = FALSE]

  list(
    matrix = expr_target,
    found = found_genes,
    missing = missing_genes
  )
}

# ============================================================================
# 3. Clinical data preparation
# ============================================================================

prepare_clinical_survival <- function(clinical) {
  cli::cli_h2("Preparing clinical survival data")

  n_initial <- nrow(clinical)

  # Derive OS.time: days_to_death if dead, otherwise days_to_last_follow_up
  clinical$OS.time <- pmax(
    as.numeric(clinical$days_to_death),
    as.numeric(clinical$days_to_last_follow_up),
    na.rm = TRUE
  )

  # Derive OS.event
  clinical$OS.event <- ifelse(
    tolower(clinical$vital_status) %in% c("dead", "deceased"),
    1, 0
  )

  # Remove missing or invalid survival times
  valid_time <- !is.na(clinical$OS.time) & clinical$OS.time > 0
  if (any(!valid_time)) {
    cli::cli_alert_info("Removing {.val {sum(!valid_time)}} patients with missing/zero/negative survival time")
    clinical <- clinical[valid_time, , drop = FALSE]
  }

  # Remove duplicated patients (keep first occurrence — already most complete from GDCquery_clinic)
  dup_ids <- duplicated(clinical$submitter_id)
  if (any(dup_ids)) {
    cli::cli_alert_info("Removing {.val {sum(dup_ids)}} duplicate clinical records")
    clinical <- clinical[!dup_ids, , drop = FALSE]
  }

  cli::cli_alert_success(
    "Clinical data: {.val {nrow(clinical)}} patients ({.val {sum(clinical$OS.event)}} events, ",
    "{.val {round(median(clinical$OS.time[clinical$OS.event == 1]) / 30.44, 1)}} months median OS)"
  )
  cli::cli_alert_info("Reduction: {.val {n_initial}} → {.val {nrow(clinical)}} patients")

  clinical
}

# ============================================================================
# 4. Expression–clinical matching
# ============================================================================

match_expression_clinical <- function(expr_matrix, clinical) {
  cli::cli_h2("Matching expression and clinical data")

  expr_patients <- colnames(expr_matrix)
  clin_patients <- clinical$submitter_id

  common <- intersect(expr_patients, clin_patients)

  cli::cli_alert_info("Expression patients: {.val {length(expr_patients)}}")
  cli::cli_alert_info("Clinical patients:   {.val {length(clin_patients)}}")
  cli::cli_alert_success("Matched patients:    {.val {length(common)}}")

  lost_expr <- setdiff(expr_patients, clin_patients)
  lost_clin <- setdiff(clin_patients, expr_patients)

  if (length(lost_expr) > 0) {
    cli::cli_alert_info("Expression-only patients (no clinical): {.val {length(lost_expr)}}")
  }
  if (length(lost_clin) > 0) {
    cli::cli_alert_info("Clinical-only patients (no expression): {.val {length(lost_clin)}}")
  }

  if (length(common) < MIN_PATIENTS_KM) {
    stop(sprintf(
      "Only %d matched patients (minimum %d required). Cannot continue survival analysis.",
      length(common), MIN_PATIENTS_KM
    ))
  }

  expr_match <- expr_matrix[, common, drop = FALSE]
  clin_match <- clinical[match(common, clinical$submitter_id), , drop = FALSE]

  list(
    expr = expr_match,
    clinical = clin_match,
    n_matched = length(common)
  )
}

# ============================================================================
# 5. Expression analysis & plots
# ============================================================================

plot_expression_boxplots <- function(expr_matrix, gene_class, outdir) {
  cli::cli_h2("Generating expression boxplots")

  genes <- rownames(expr_matrix)
  n_patients <- ncol(expr_matrix)

  # Build long-format data
  expr_long <- data.frame(
    gene = factor(rep(genes, each = n_patients), levels = genes),
    sample = rep(colnames(expr_matrix), times = length(genes)),
    expression = as.vector(t(expr_matrix)),
    class = rep(gene_class[genes], each = n_patients),
    stringsAsFactors = FALSE
  )

  # Add class-aware ordering
  expr_long$gene <- factor(expr_long$gene, levels = genes)

  # Color by gene class
  class_colors <- gene_class_palette[gene_class[genes]]

  p <- ggplot(expr_long, aes(x = gene, y = expression, fill = class)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.5, width = 0.6,
                 position = position_dodge(0.8)) +
    scale_fill_manual(values = gene_class_palette, name = "Gene class") +
    labs(
      title = "TCGA-GBM: Target Gene Expression",
      subtitle = sprintf("Primary tumor samples (n = %d patients)", n_patients),
      x = NULL,
      y = expression(log[2] * "(TPM + 1)")
    ) +
    theme_gbm_journal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "right"
    )

  out <- file.path(outdir, "Fig01A_expression_boxplot")
  save_figure_pub(p, out, width = 10, height = 6)
  p
}

plot_expression_violin <- function(expr_matrix, gene_class, outdir) {
  cli::cli_h2("Generating faceted violin plots")

  genes <- rownames(expr_matrix)
  n_patients <- ncol(expr_matrix)

  expr_long <- data.frame(
    gene = rep(genes, each = n_patients),
    sample = rep(colnames(expr_matrix), times = length(genes)),
    expression = as.vector(t(expr_matrix)),
    class = rep(gene_class[genes], each = n_patients),
    stringsAsFactors = FALSE
  )

  expr_long$gene <- factor(expr_long$gene, levels = genes)

  p <- ggplot(expr_long, aes(x = gene, y = expression, fill = class)) +
    geom_violin(trim = TRUE, alpha = 0.7, width = 0.8) +
    geom_boxplot(width = 0.15, outlier.size = 0.4, alpha = 0.8,
                 position = position_dodge(0.8)) +
    geom_jitter(width = 0.08, size = 0.3, alpha = 0.25, color = "#333333") +
    scale_fill_manual(values = gene_class_palette, name = "Gene class") +
    facet_wrap(~ gene, ncol = 4, scales = "free_y") +
    labs(
      title = "TCGA-GBM: Gene Expression Distribution",
      subtitle = sprintf("Primary tumors, n = %d patients. Violin + boxplot + jitter", n_patients),
      x = NULL,
      y = expression(log[2] * "(TPM + 1)")
    ) +
    theme_gbm_journal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(face = "bold", size = 9)
    )

  out <- file.path(outdir, "Fig01A2_expression_violin_faceted")
  save_figure_pub(p, out, width = 12, height = 10)
  p
}

plot_expression_heatmap <- function(expr_matrix, gene_class, outdir) {
  n_genes <- nrow(expr_matrix)

  if (n_genes < MIN_GENES_HEATMAP) {
    cli::cli_alert_info(
      "Skipping heatmap: {.val {n_genes}} genes found (< {.val {MIN_GENES_HEATMAP}} threshold). ",
      "Heatmap not meaningful with so few genes."
    )
    return(invisible(NULL))
  }

  cli::cli_h2("Generating expression heatmap")

  # Z-score normalize per gene
  expr_z <- t(scale(t(expr_matrix)))

  # Annotation
  row_anno <- ComplexHeatmap::rowAnnotation(
    `Gene class` = gene_class[rownames(expr_z)],
    col = list(`Gene class` = gene_class_palette),
    annotation_legend_param = list(
      `Gene class` = list(title = "Gene class")
    )
  )

  # Color mapping
  col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c(viridis::viridis(3, direction = -1))
  )

  ht <- ComplexHeatmap::Heatmap(
    expr_z,
    name = "Z-score",
    col = col_fun,
    left_annotation = row_anno,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    column_title = sprintf("TCGA-GBM Patients (n = %d)", ncol(expr_z)),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 9),
    heatmap_legend_param = list(direction = "horizontal")
  )

  out <- file.path(outdir, "Fig01B_expression_heatmap")
  pdf(paste0(out, ".pdf"), width = 14, height = 8)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  dev.off()
  png(paste0(out, ".png"), width = 14, height = 8, units = "in", res = 600)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  dev.off()
  tiff(paste0(out, ".tiff"), width = 14, height = 8, units = "in", res = 600, compression = "lzw")
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  dev.off()

  cli::cli_alert_success("{.file Fig01B_expression_heatmap} → .pdf, .png, .tiff")
  invisible(ht)
}

build_expression_summary_table <- function(expr_matrix, gene_class, found_genes, missing_genes) {
  cli::cli_h2("Building expression summary table")

  genes <- rownames(expr_matrix)
  n_patients <- ncol(expr_matrix)

  df <- data.frame(
    gene = genes,
    gene_class = gene_class[genes],
    found_in_tcga = TRUE,
    mean_expression = round(rowMeans(expr_matrix, na.rm = TRUE), 3),
    median_expression = round(apply(expr_matrix, 1, median, na.rm = TRUE), 3),
    sd_expression = round(apply(expr_matrix, 1, sd, na.rm = TRUE), 3),
    min_expression = round(apply(expr_matrix, 1, min, na.rm = TRUE), 3),
    max_expression = round(apply(expr_matrix, 1, max, na.rm = TRUE), 3),
    n_samples = n_patients,
    n_patients = n_patients,
    stringsAsFactors = FALSE
  )

  # Add missing genes
  if (length(missing_genes) > 0) {
    df_missing <- data.frame(
      gene = missing_genes,
      gene_class = gene_class[missing_genes],
      found_in_tcga = FALSE,
      mean_expression = NA_real_,
      median_expression = NA_real_,
      sd_expression = NA_real_,
      min_expression = NA_real_,
      max_expression = NA_real_,
      n_samples = 0L,
      n_patients = 0L,
      stringsAsFactors = FALSE
    )
    df <- rbind(df, df_missing)
    df <- df[match(c(genes, missing_genes), df$gene), , drop = FALSE]
  }

  # Sort by mean expression (descending), NA last
  df <- df[order(-df$mean_expression, na.last = TRUE), , drop = FALSE]
  rownames(df) <- NULL

  df
}

# ============================================================================
# 6. Survival analysis
# ============================================================================

run_km_for_gene <- function(gene, expr_vec, clinical, min_events = MIN_EVENTS_PER_GROUP) {
  # Median split
  median_expr <- median(expr_vec, na.rm = TRUE)
  group <- ifelse(expr_vec > median_expr, "High expression", "Low expression")
  group <- factor(group, levels = c("Low expression", "High expression"))

  n_high <- sum(group == "High expression")
  n_low  <- sum(group == "Low expression")

  # Per-gene survival data frame. Using a column-based formula with `data=`
  # (instead of a bare local Surv object) keeps the fitted model objects
  # self-contained, so survminer can re-evaluate the formula downstream
  # without hitting "object 'surv_obj' not found".
  surv_df <- data.frame(
    OS.time = clinical$OS.time,
    OS.event = clinical$OS.event,
    group = group,
    expr = expr_vec
  )

  # Log-rank test
  surv_diff <- survdiff(Surv(OS.time, OS.event) ~ group, data = surv_df)
  logrank_p <- 1 - pchisq(surv_diff$chisq, df = 1)

  # Kaplan-Meier fit
  surv_fit <- survfit(Surv(OS.time, OS.event) ~ group, data = surv_df)

  # Cox proportional hazards
  cox_fit <- coxph(Surv(OS.time, OS.event) ~ group, data = surv_df)
  cox_summary <- summary(cox_fit)
  cox_hr <- cox_summary$conf.int[1, "exp(coef)"]
  cox_ci_lower <- cox_summary$conf.int[1, "lower .95"]
  cox_ci_upper <- cox_summary$conf.int[1, "upper .95"]
  cox_p <- cox_summary$coefficients[1, "Pr(>|z|)"]

  # PH assumption check
  zph_test <- tryCatch(
    cox.zph(cox_fit),
    error = function(e) NULL
  )
  ph_p <- if (!is.null(zph_test)) zph_test$table[1, "p"] else NA_real_

  # Median OS per group
  median_os_high <- if (any(group == "High expression")) {
    fit_high <- survfit(Surv(OS.time, OS.event) ~ 1,
                        data = surv_df[surv_df$group == "High expression", ])
    summary(fit_high)$table["median"] / 30.44  # convert days to months
  } else NA_real_
  median_os_low <- if (any(group == "Low expression")) {
    fit_low <- survfit(Surv(OS.time, OS.event) ~ 1,
                       data = surv_df[surv_df$group == "Low expression", ])
    summary(fit_low)$table["median"] / 30.44
  } else NA_real_

  n_events <- sum(clinical$OS.event)

  list(
    gene = gene,
    cutoff_method = "Median",
    cutoff_value = median_expr,
    n_patients = length(expr_vec),
    n_events = n_events,
    high_group_n = n_high,
    low_group_n = n_low,
    median_OS_high_months = round(median_os_high, 1),
    median_OS_low_months = round(median_os_low, 1),
    logrank_p = logrank_p,
    cox_HR = cox_hr,
    cox_CI_lower = cox_ci_lower,
    cox_CI_upper = cox_ci_upper,
    cox_p = cox_p,
    PH_assumption_p = ph_p,
    surv_diff = surv_diff,
    surv_fit = surv_fit,
    cox_fit = cox_fit,
    group = group,
    df = surv_df,
    clin = clinical
  )
}

plot_km_gene <- function(km_result, gene_class, outdir) {
  gene <- km_result$gene
  class_label <- gene_class[gene]
  if (is.na(class_label)) class_label <- "Unknown"

  hr_str <- sprintf("%.2f", km_result$cox_HR)
  ci_str <- sprintf("%.2f–%.2f", km_result$cox_CI_lower, km_result$cox_CI_upper)
  p_str <- format_pval(km_result$logrank_p)

  # Build survfit for ggsurvplot from the gene's self-contained data frame
  fit <- survfit(Surv(OS.time, OS.event) ~ group, data = km_result$df)

  p <- ggsurvplot(
    fit,
    data = km_result$df,
    pval = FALSE,  # we add our own
    risk.table = TRUE,
    risk.table.height = 0.25,
    censor.shape = "+",
    censor.size = 2,
    palette = km_group_palette,
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    title = sprintf("%s (%s) — TCGA-GBM", gene, class_label),
    subtitle = sprintf("HR = %s (95%% CI %s), log-rank P = %s | n = %d, events = %d",
                       hr_str, ci_str, p_str,
                       km_result$n_patients, km_result$n_events),
    legend.title = "",
    legend.labs = c("Low expression", "High expression"),
    ggtheme = theme_gbm_journal(base_size = 10),
    risk.table.y.text = TRUE,
    surv.median.line = "hv",
    xlim = c(0, max(km_result$clin$OS.time, na.rm = TRUE) / 30.44 * 1.05),
    break.x.by = 12  # tick every 12 months
  )

  # Convert x-axis days → months in the plot object
  p$plot <- p$plot +
    scale_x_continuous(
      breaks = seq(0, ceiling(max(km_result$clin$OS.time) / 30.44), by = 12),
      labels = seq(0, ceiling(max(km_result$clin$OS.time) / 30.44), by = 12)
    )

  # Also fix risk table x-axis
  if (!is.null(p$table)) {
    p$table <- p$table +
      scale_x_continuous(
        breaks = seq(0, ceiling(max(km_result$clin$OS.time) / 30.44), by = 12),
        labels = seq(0, ceiling(max(km_result$clin$OS.time) / 30.44), by = 12)
      )
  }

  out <- file.path(outdir, sprintf("Fig01C_km_%s", gene))
  save_figure_km(p, out, width = 9, height = 8)

  invisible(p)
}

plot_km_combined <- function(km_results, gene_class, outdir, ncol = 4) {
  cli::cli_h2("Generating combined KM figure")

  n_genes <- length(km_results)
  if (n_genes == 0) {
    cli::cli_alert_warning("No KM results to plot")
    return(invisible(NULL))
  }

  # Extract individual ggplot objects from ggsurvplot results
  plot_list <- lapply(names(km_results), function(gene) {
    res <- km_results[[gene]]
    class_label <- gene_class[gene]
    if (is.na(class_label)) class_label <- "Unknown"

    fit <- survfit(Surv(OS.time, OS.event) ~ group, data = res$df)

    hr_str <- sprintf("%.2f", res$cox_HR)
    p_str <- format_pval(res$logrank_p)

    # Minimal KM plot for combined figure
    ggp <- ggsurvplot(
      fit,
      data = res$df,
      pval = FALSE,
      risk.table = FALSE,
      censor.shape = "+",
      censor.size = 1.5,
      palette = km_group_palette,
      xlab = NULL,
      ylab = NULL,
      title = sprintf("%s (%s)", gene, class_label),
      legend = "none",
      ggtheme = theme_gbm_journal(base_size = 8),
      surv.median.line = "hv",
      xlim = c(0, max(res$clin$OS.time, na.rm = TRUE) / 30.44 * 1.02)
    )

    # Add annotation
    ggp$plot +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("P = %s", p_str),
               hjust = 1.1, vjust = 1.5, size = 2.8) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("HR = %s", hr_str),
               hjust = 1.1, vjust = 3.2, size = 2.8)
  })

  # Arrange with patchwork
  n_rows <- ceiling(n_genes / ncol)
  combined <- patchwork::wrap_plots(plot_list, ncol = ncol, nrow = n_rows) +
    patchwork::plot_annotation(
      title = "TCGA-GBM: Overall Survival by Gene Expression",
      caption = sprintf("Median-split high vs low expression. n = %d patients per gene.",
                        km_results[[1]]$n_patients)
    )

  out <- file.path(outdir, "Fig01C_km_combined")
  pdf(paste0(out, ".pdf"), width = 16, height = 4 * n_rows)
  print(combined)
  dev.off()
  png(paste0(out, ".png"), width = 16, height = 4 * n_rows, units = "in", res = 600)
  print(combined)
  dev.off()
  tiff(paste0(out, ".tiff"), width = 16, height = 4 * n_rows, units = "in", res = 600, compression = "lzw")
  print(combined)
  dev.off()

  cli::cli_alert_success("{.file Fig01C_km_combined} → .pdf, .png, .tiff")
  invisible(combined)
}

plot_survival_forest <- function(survival_summary, outdir) {
  cli::cli_h2("Generating forest plot")

  df <- survival_summary[order(survival_summary$cox_HR), , drop = FALSE]
  df$gene <- factor(df$gene, levels = df$gene)

  # Determine significance coloring
  df$sig <- ifelse(df$cox_p < 0.05,
                   ifelse(df$cox_HR < 1, "Better survival (P < 0.05)", "Worse survival (P < 0.05)"),
                   "Not significant")

  p <- ggplot(df, aes(x = cox_HR, y = gene, color = gene_class)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#999999", linewidth = 0.5) +
    geom_point(aes(size = n_patients), alpha = 0.9) +
    geom_errorbarh(aes(xmin = cox_CI_lower, xmax = cox_CI_upper),
                   height = 0.2, linewidth = 0.7) +
    scale_x_log10(
      breaks = c(0.5, 0.7, 1.0, 1.4, 2.0),
      labels = c("0.5", "0.7", "1.0", "1.4", "2.0"),
      limits = c(
        max(0.3, min(df$cox_CI_lower, na.rm = TRUE) * 0.8),
        min(3.0, max(df$cox_CI_upper, na.rm = TRUE) * 1.2)
      )
    ) +
    scale_color_manual(values = gene_class_palette, name = "Gene class") +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    labs(
      title = "TCGA-GBM: Cox Regression Forest Plot",
      subtitle = sprintf("Overall survival, median-split expression. n = %d patients per gene.",
                         df$n_patients[1]),
      x = "Hazard Ratio (log scale)",
      y = NULL
    ) +
    theme_gbm_journal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    ) +
    geom_text(
      aes(label = sprintf("p=%.3f", cox_p), x = cox_CI_upper),
      hjust = -0.3, size = 3, color = "#444444"
    )

  out <- file.path(outdir, "Fig01D_survival_forestplot")
  save_figure_pub(p, out, width = 10, height = 6)
  p
}

# ============================================================================
# 6b. Advanced / robustness survival analyses (Phase 2)
#
# The PRIMARY, pre-specified analysis is the median split (Table02). The
# functions below are SUPPORTING / robustness analyses:
#   - Continuous Cox  : uses the full expression range (no dichotomization)
#   - Multivariable Cox: median-split group adjusted for age + sex
#   - Optimal cutpoint : EXPLORATORY only (optimism-biased; inflates p-values)
# Primary and exploratory results are kept in separate columns and never
# conflated. Nulls are reported as nulls.
# ============================================================================

run_advanced_survival <- function(km_results, clinical, surv_summary) {
  cli::cli_h2("Advanced survival analyses (continuous / multivariable / exploratory cutpoint)")

  # Covariate availability. TCGA-GBM GDC clinical exposes age + sex only;
  # IDH mutation and MGMT methylation are NOT in this table.
  has_age <- "age_at_index" %in% names(clinical)
  has_sex <- "gender" %in% names(clinical)
  covar_terms <- c(if (has_age) "age", if (has_sex) "sex")
  covar_label <- if (length(covar_terms)) paste(covar_terms, collapse = " + ") else "none"
  cli::cli_alert_info("Multivariable covariates: {.val {covar_label}} (IDH/MGMT not in GDC clinical)")

  rows <- lapply(names(km_results), function(gene) {
    r <- km_results[[gene]]
    d <- r$df  # columns: OS.time, OS.event, group, expr (rows aligned to `clinical`)

    # --- Continuous Cox (HR per log2 unit and per SD) ---
    cont <- tryCatch({
      cc <- coxph(Surv(OS.time, OS.event) ~ expr, data = d)
      sc <- summary(cc)
      b  <- sc$coefficients[1, "coef"]
      se <- sc$coefficients[1, "se(coef)"]
      sdv <- sd(d$expr, na.rm = TRUE)
      list(hr_log2 = exp(b),
           hr_sd = exp(b * sdv),
           ci_sd_l = exp((b - 1.96 * se) * sdv),
           ci_sd_u = exp((b + 1.96 * se) * sdv),
           p = sc$coefficients[1, "Pr(>|z|)"])
    }, error = function(e) NULL)

    # --- Multivariable Cox: group + age + sex ---
    mv <- tryCatch({
      md <- d
      if (has_age) md$age <- suppressWarnings(as.numeric(clinical$age_at_index))
      if (has_sex) md$sex <- factor(clinical$gender)
      rhs <- paste(c("group", covar_terms), collapse = " + ")
      mfit <- coxph(as.formula(paste("Surv(OS.time, OS.event) ~", rhs)), data = md)
      ms <- summary(mfit)
      grow <- grep("^group", rownames(ms$coefficients))[1]
      list(hr = ms$conf.int[grow, "exp(coef)"],
           ci_l = ms$conf.int[grow, "lower .95"],
           ci_u = ms$conf.int[grow, "upper .95"],
           p = ms$coefficients[grow, "Pr(>|z|)"],
           n = mfit$n)
    }, error = function(e) NULL)

    # --- Exploratory optimal cutpoint (survminer::surv_cutpoint) ---
    cut <- tryCatch({
      cdf <- data.frame(OS.time = d$OS.time, OS.event = d$OS.event, expr = d$expr)
      scut <- survminer::surv_cutpoint(cdf, time = "OS.time", event = "OS.event",
                                       variables = "expr", minprop = 0.1)
      cutval <- scut$cutpoint$cutpoint[1]
      cat_df <- survminer::surv_categorize(scut)
      sd2 <- survdiff(Surv(OS.time, OS.event) ~ expr, data = cat_df)
      list(value = cutval,
           high_n = sum(cat_df$expr == "high"),
           low_n  = sum(cat_df$expr == "low"),
           p = 1 - pchisq(sd2$chisq, 1))
    }, error = function(e) NULL)

    data.frame(
      gene = gene,
      gene_class = GENE_CLASS[gene],
      primary_logrank_p = surv_summary$logrank_p[match(gene, surv_summary$gene)],
      primary_adj_p     = surv_summary$adjusted_p[match(gene, surv_summary$gene)],
      cont_HR_per_log2  = if (!is.null(cont)) round(cont$hr_log2, 3) else NA_real_,
      cont_HR_per_SD    = if (!is.null(cont)) round(cont$hr_sd, 3) else NA_real_,
      cont_CI_SD_lower  = if (!is.null(cont)) round(cont$ci_sd_l, 3) else NA_real_,
      cont_CI_SD_upper  = if (!is.null(cont)) round(cont$ci_sd_u, 3) else NA_real_,
      cont_p            = if (!is.null(cont)) cont$p else NA_real_,
      mv_HR_group       = if (!is.null(mv)) round(mv$hr, 3) else NA_real_,
      mv_CI_lower       = if (!is.null(mv)) round(mv$ci_l, 3) else NA_real_,
      mv_CI_upper       = if (!is.null(mv)) round(mv$ci_u, 3) else NA_real_,
      mv_p              = if (!is.null(mv)) mv$p else NA_real_,
      mv_n              = if (!is.null(mv)) mv$n else NA_integer_,
      mv_covariates     = covar_label,
      cutpoint_value    = if (!is.null(cut)) round(cut$value, 3) else NA_real_,
      cutpoint_high_n   = if (!is.null(cut)) cut$high_n else NA_integer_,
      cutpoint_low_n    = if (!is.null(cut)) cut$low_n else NA_integer_,
      cutpoint_logrank_p_EXPLORATORY = if (!is.null(cut)) cut$p else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  adv <- do.call(rbind, rows)
  adv$cont_adj_p <- p.adjust(adv$cont_p, method = "BH")  # FDR for continuous Cox

  ord <- c("gene", "gene_class", "primary_logrank_p", "primary_adj_p",
           "cont_HR_per_log2", "cont_HR_per_SD", "cont_CI_SD_lower", "cont_CI_SD_upper",
           "cont_p", "cont_adj_p",
           "mv_HR_group", "mv_CI_lower", "mv_CI_upper", "mv_p", "mv_n", "mv_covariates",
           "cutpoint_value", "cutpoint_high_n", "cutpoint_low_n",
           "cutpoint_logrank_p_EXPLORATORY")
  adv <- adv[order(adv$cont_p), ord]
  adv
}

plot_continuous_cox_forest <- function(adv_df, outdir) {
  cli::cli_h2("Generating continuous-Cox forest plot")
  df <- adv_df[!is.na(adv_df$cont_HR_per_SD), , drop = FALSE]
  if (nrow(df) == 0) {
    cli::cli_alert_warning("No continuous Cox results to plot")
    return(invisible(NULL))
  }
  df <- df[order(df$cont_HR_per_SD), ]
  df$gene <- factor(df$gene, levels = df$gene)

  p <- ggplot(df, aes(x = cont_HR_per_SD, y = gene, color = gene_class)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#999999", linewidth = 0.5) +
    geom_point(size = 3, alpha = 0.9) +
    geom_errorbarh(aes(xmin = cont_CI_SD_lower, xmax = cont_CI_SD_upper),
                   height = 0.2, linewidth = 0.7) +
    scale_x_log10() +
    scale_color_manual(values = gene_class_palette, name = "Gene class") +
    labs(
      title = "TCGA-GBM: Continuous Cox (HR per 1 SD of log2 expression)",
      subtitle = "Robustness analysis using the full expression range (no dichotomization)",
      x = "Hazard Ratio per SD (log scale)",
      y = NULL
    ) +
    theme_gbm_journal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "right") +
    geom_text(aes(label = sprintf("p=%.3f", cont_p), x = cont_CI_SD_upper),
              hjust = -0.3, size = 3, color = "#444444")

  out <- file.path(outdir, "Fig01E_continuous_cox_forest")
  save_figure_pub(p, out, width = 10, height = 6)
  invisible(p)
}

write_robustness_notes <- function(adv_df, outdir) {
  n_g       <- nrow(adv_df)
  n_cont    <- sum(adv_df$cont_p < 0.05, na.rm = TRUE)
  n_cont_fdr <- sum(adv_df$cont_adj_p < 0.05, na.rm = TRUE)
  n_mv      <- sum(adv_df$mv_p < 0.05, na.rm = TRUE)
  n_cut     <- sum(adv_df$cutpoint_logrank_p_EXPLORATORY < 0.05, na.rm = TRUE)

  lines <- c(
    "# Module 01 — Survival Robustness Analyses (Phase 2)",
    "",
    sprintf("**Generated:** %s", format(Sys.time())),
    "",
    "## Purpose",
    "Supporting analyses that test the robustness of the PRIMARY median-split",
    "overall-survival results (Table02). They do NOT replace the primary analysis.",
    "",
    "## Methods",
    "- **Continuous Cox** — HR per 1 SD of log2(TPM+1) over the full expression",
    "  range (no dichotomization). BH-FDR applied across genes (`cont_adj_p`).",
    "- **Multivariable Cox** — median-split group adjusted for age + sex.",
    "  IDH mutation, MGMT methylation, and molecular subtype are NOT available in",
    "  the GDC clinical table for TCGA-GBM and could not be included.",
    "- **Optimal cutpoint (EXPLORATORY)** — `survminer::surv_cutpoint` selects the",
    "  split that maximizes the log-rank statistic. This is optimism-biased: the",
    "  resulting p-values (`cutpoint_logrank_p_EXPLORATORY`) are NOT corrected for",
    "  cutpoint selection and are hypothesis-generating ONLY.",
    "",
    "## Results summary",
    sprintf("- Continuous Cox: %d/%d genes nominally significant (p<0.05); %d survive FDR.",
            n_cont, n_g, n_cont_fdr),
    sprintf("- Multivariable Cox (group | age + sex): %d/%d genes with significant group term.",
            n_mv, n_g),
    sprintf("- Exploratory optimal cutpoint: %d/%d genes nominally significant (EXPLORATORY, uncorrected).",
            n_cut, n_g),
    "",
    "## Interpretation caveats",
    "- 'Essential in a CRISPR depletion screen' (A172-S vs A172-R) is a different",
    "  biological claim from 'prognostic for overall survival in TCGA-GBM patients'.",
    "  Concordance between the two is not expected a priori.",
    "- Cutpoint-optimized p-values must never be reported as primary evidence.",
    "- See Table04_survival_robustness for the full per-gene values."
  )
  writeLines(lines, file.path(outdir, "01_survival_robustness_notes.md"))
  cli::cli_alert_success("01_survival_robustness_notes.md written")
}

build_survival_summary_table <- function(km_results) {
  cli::cli_h2("Building survival summary table")

  df <- do.call(rbind, lapply(km_results, function(r) {
    data.frame(
      gene = r$gene,
      cutoff_method = r$cutoff_method,
      cutoff_value = round(r$cutoff_value, 3),
      n_patients = r$n_patients,
      n_events = r$n_events,
      high_group_n = r$high_group_n,
      low_group_n = r$low_group_n,
      median_OS_high_months = r$median_OS_high_months,
      median_OS_low_months = r$median_OS_low_months,
      logrank_p = r$logrank_p,
      cox_HR = round(r$cox_HR, 3),
      cox_CI_lower = round(r$cox_CI_lower, 3),
      cox_CI_upper = round(r$cox_CI_upper, 3),
      cox_p = r$cox_p,
      PH_assumption_p = r$PH_assumption_p,
      stringsAsFactors = FALSE
    )
  }))

  # Add gene class
  df$gene_class <- GENE_CLASS[df$gene]

  # Adjusted p-value (Benjamini-Hochberg FDR)
  df$adjusted_p <- p.adjust(df$logrank_p, method = "BH")

  # Interpretation
  df$interpretation <- NA_character_
  for (i in seq_len(nrow(df))) {
    if (df$cox_p[i] < 0.05 && df$cox_HR[i] < 1) {
      df$interpretation[i] <- "Significant association with better survival"
    } else if (df$cox_p[i] < 0.05 && df$cox_HR[i] > 1) {
      df$interpretation[i] <- "Significant association with worse survival"
    } else if (df$cox_p[i] < 0.10) {
      df$interpretation[i] <- "Near-significant association (exploratory)"
    } else {
      df$interpretation[i] <- "No significant association"
    }
  }

  # Flag PH assumption violations
  ph_violated <- which(df$PH_assumption_p < 0.05)
  if (length(ph_violated) > 0) {
    df$interpretation[ph_violated] <- paste(
      df$interpretation[ph_violated],
      "[PH assumption violated — interpret with caution]"
    )
  }

  # Reorder columns
  df <- df[, c("gene", "gene_class", "cutoff_method", "cutoff_value",
               "n_patients", "n_events", "high_group_n", "low_group_n",
               "median_OS_high_months", "median_OS_low_months",
               "logrank_p", "cox_HR", "cox_CI_lower", "cox_CI_upper",
               "cox_p", "adjusted_p", "PH_assumption_p", "interpretation")]

  # Sort by log-rank p
  df <- df[order(df$logrank_p), , drop = FALSE]
  rownames(df) <- NULL

  df
}

# ============================================================================
# 7. QC table & summary
# ============================================================================

build_sample_qc_table <- function(n_expression_aliquots, n_unique_patients,
                                   n_tumor_samples, n_normal_samples,
                                   n_matched_clinical, found_genes, missing_genes,
                                   excluded_details = NULL) {
  cli::cli_h2("Building sample QC table")

  n_total_genes <- length(found_genes) + length(missing_genes)

  exclusion_reasons <- c()
  if (length(missing_genes) > 0) {
    exclusion_reasons <- c(
      exclusion_reasons,
      sprintf("%d genes not found in TCGA: %s", length(missing_genes),
              paste(missing_genes, collapse = ", "))
    )
  }
  if (n_normal_samples < MIN_NORMAL_SAMPLES) {
    exclusion_reasons <- c(
      exclusion_reasons,
      sprintf("Tumor-normal comparison skipped (%d normal samples < %d minimum)",
              n_normal_samples, MIN_NORMAL_SAMPLES)
    )
  }

  df <- data.frame(
    metric = c(
      "Total expression samples (aliquots)",
      "Unique patients (after averaging)",
      "Primary tumor samples",
      "Solid tissue normal samples",
      "Target genes queried",
      "Genes found in TCGA",
      "Genes missing from TCGA",
      "Matched expression + clinical patients",
      "Excluded samples/genes",
      "Tumor-normal comparison",
      "KM split method",
      "Survival endpoint"
    ),
    value = c(
      as.character(n_expression_aliquots),
      as.character(n_unique_patients),
      as.character(n_tumor_samples),
      as.character(n_normal_samples),
      as.character(n_total_genes),
      as.character(length(found_genes)),
      as.character(length(missing_genes)),
      as.character(n_matched_clinical),
      if (length(exclusion_reasons) > 0) paste(exclusion_reasons, collapse = "; ") else "None",
      if (n_normal_samples >= MIN_NORMAL_SAMPLES) "Tumor vs Normal" else "Tumor only (insufficient normals)",
      "Median expression split (High vs Low)",
      "Overall Survival (OS)"
    ),
    stringsAsFactors = FALSE
  )

  df
}

write_module_summary <- function(expr_summary, surv_summary, sample_qc,
                                  found_genes, missing_genes,
                                  n_matched, n_normal, outdir) {
  cli::cli_h2("Writing module summary")

  # Count significant genes
  sig_genes <- surv_summary[surv_summary$cox_p < 0.05, ]
  near_sig <- surv_summary[surv_summary$cox_p >= 0.05 & surv_summary$cox_p < 0.10, ]
  hr_gt_1 <- surv_summary[surv_summary$cox_HR > 1, ]
  hr_lt_1 <- surv_summary[surv_summary$cox_HR < 1, ]

  # Build summary
  md_lines <- c(
    sprintf("# Module 01: TCGA-GBM Expression & Survival Analysis"),
    sprintf(""),
    sprintf("**Generated:** %s", Sys.time()),
    sprintf("**TCGA Project:** TCGA-GBM (Glioblastoma Multiforme)"),
    sprintf("**Data:** RNA-seq (STAR-Counts), Clinical (GDC)"),
    sprintf(""),
    sprintf("## 1. Data Acquisition Summary"),
    sprintf(""),
    sprintf("- **Expression samples:** %d aliquots → %d unique tumor patients",
            nrow(sample_qc) > 0, ncol(sample_qc)),
    sprintf("- **Normal samples:** %d (tumor-normal comparison %s)",
            n_normal,
            if (n_normal >= MIN_NORMAL_SAMPLES) "performed" else "skipped — insufficient samples"),
    sprintf("- **Matched expression + clinical:** %d patients", n_matched),
    sprintf(""),
    sprintf("## 2. Gene Detection"),
    sprintf(""),
    sprintf("### Found (%d/%d genes):", length(found_genes), length(found_genes) + length(missing_genes)),
    sprintf(""),
    paste("-", found_genes, collapse = "\n"),
    sprintf(""),
    if (length(missing_genes) > 0) {
      c(
        sprintf("### Missing (%d genes):", length(missing_genes)),
        sprintf(""),
        paste("-", missing_genes, collapse = "\n"),
        sprintf(""),
        sprintf("**Reason:** Gene symbol not present in TCGA-GBM STAR-Counts gene annotation."),
        sprintf("")
      )
    } else {
      c(sprintf("All target genes were found in TCGA-GBM expression data."), "")
    },
    sprintf("## 3. Expression Patterns"),
    sprintf(""),
    sprintf("| Gene | Class | Mean log2(TPM+1) | Median log2(TPM+1) |"),
    sprintf("|------|-------|-------------------|--------------------|"),
    paste(sprintf("| %s | %s | %.2f | %.2f |",
                  expr_summary$gene, expr_summary$gene_class,
                  expr_summary$mean_expression, expr_summary$median_expression),
          collapse = "\n"),
    sprintf(""),
    sprintf("## 4. Survival Analysis"),
    sprintf(""),
    sprintf("**Method:** Median expression split, overall survival, Cox proportional hazards."),
    sprintf(""),
    if (nrow(sig_genes) > 0) {
      c(
        sprintf("### Genes with Significant Survival Association (P < 0.05):"),
        sprintf(""),
        paste(sprintf("- **%s** (%s): HR = %.2f [%.2f–%.2f], P = %.3f",
                      sig_genes$gene, sig_genes$gene_class,
                      sig_genes$cox_HR, sig_genes$cox_CI_lower, sig_genes$cox_CI_upper,
                      sig_genes$cox_p),
              collapse = "\n"),
        sprintf("")
      )
    } else {
      c(sprintf("No genes showed statistically significant survival association (P < 0.05)."), "")
    },
    if (nrow(near_sig) > 0) {
      c(
        sprintf("### Genes with Near-Significant Association (P < 0.10) [EXPLORATORY]:"),
        sprintf(""),
        paste(sprintf("- **%s** (%s): HR = %.2f, P = %.3f",
                      near_sig$gene, near_sig$gene_class,
                      near_sig$cox_HR, near_sig$cox_p),
              collapse = "\n"),
        sprintf("")
      )
    },
    sprintf("### Direction of Effect:"),
    sprintf("- **HR > 1 (high expression = worse survival):** %d genes", nrow(hr_gt_1)),
    paste(" ", "-", sprintf("%s (HR=%.2f)", hr_gt_1$gene, hr_gt_1$cox_HR), collapse = "\n"),
    sprintf("- **HR < 1 (high expression = better survival):** %d genes", nrow(hr_lt_1)),
    paste(" ", "-", sprintf("%s (HR=%.2f)", hr_lt_1$gene, hr_lt_1$cox_HR), collapse = "\n"),
    sprintf(""),
    sprintf("## 5. CRISPR Screen Validation Assessment"),
    sprintf(""),
    sprintf("### Caveats"),
    sprintf(""),
    sprintf("1. **This is external clinical validation** — TCGA-GBM data reflect tumor biology"),
    sprintf("   in treatment-naïve patients, not CRISPR-engineered cell lines."),
    sprintf("2. **The CRISPR screen identified genes whose depletion confers TRAIL resistance**"),
    sprintf("   in A172-S vs A172-R cells. Survival associations in TCGA-GBM may reflect"),
    sprintf("   different biological mechanisms."),
    sprintf("3. **Sample size** — TCGA-GBM typically has ~150 patients with expression data,"),
    sprintf("   providing limited power for survival analysis."),
    sprintf("4. **Cox models are univariate** — no adjustment for age, sex, IDH status, or MGMT methylation."),
    sprintf("5. **Exploratory results** — P-values are unadjusted for multiple testing across genes."),
    sprintf("   Adjusted (FDR) p-values are provided in Table02."),
    sprintf(""),
    sprintf("### Summary"),
    sprintf(""),
    sprintf("All findings should be interpreted as hypothesis-generating. "),
    sprintf("Significant survival associations in TCGA-GBM do not directly validate "),
    sprintf("the CRISPR screen mechanism. They suggest that these metabolic genes "),
    sprintf("may have prognostic relevance in glioblastoma, warranting further "),
    sprintf("mechanistic investigation."),
    sprintf(""),
    sprintf("---"),
    sprintf(""),
    sprintf("*This summary was automatically generated by 01_expression_survival.R v2.0.*"),
    sprintf("*Do not overclaim statistical significance. All results are exploratory.*")
  )

  out <- file.path(outdir, "01_expression_survival_summary.md")
  writeLines(md_lines, out)
  cli::cli_alert_success("{.file 01_expression_survival_summary.md} written")

  invisible(md_lines)
}

export_session_info <- function(outdir) {
  cli::cli_h2("Exporting session info")
  si <- capture.output(sessionInfo())
  out <- file.path(outdir, "sessionInfo.txt")
  writeLines(si, out)
  cli::cli_alert_success("{.file sessionInfo.txt} written")
  invisible(si)
}

# ============================================================================
# 8. Main pipeline
# ============================================================================

run_pipeline <- function() {
  start_time <- Sys.time()

  # --- 0. Setup ---
  load_packages()
  cli::cli_h1("Pipeline Execution")

  # Create output directories
  dir.create(FIGURES_DIR_01, showWarnings = FALSE, recursive = TRUE)
  dir.create(TABLES_DIR_01, showWarnings = FALSE, recursive = TRUE)
  cli::cli_alert_info("Output: {.path {OUTPUT_DIR_01}}")

  # --- 1. Download data ---
  se <- download_tcga_expression(TCGA_PROJECT, CACHE_DIR)
  clinical_raw <- download_tcga_clinical(TCGA_PROJECT)

  # --- 2. Prepare expression ---
  expr_all <- prepare_expression_matrix(se)
  sample_info <- detect_sample_types(se)

  # Select tumor samples and average by patient
  expr_tumor_patient <- select_and_average_tumor_samples(expr_all, sample_info)

  # Extract target genes
  target_result <- extract_target_genes(expr_tumor_patient, ALL_GENES)
  expr_target <- target_result$matrix
  found_genes <- target_result$found
  missing_genes <- target_result$missing

  if (length(found_genes) == 0) {
    stop("No target genes found in TCGA-GBM. Cannot continue.")
  }

  # --- 3. Prepare clinical ---
  clinical <- prepare_clinical_survival(clinical_raw)

  # --- 4. Match ---
  matched <- match_expression_clinical(expr_target, clinical)

  # --- 5. Expression figures ---
  boxplot_p <- plot_expression_boxplots(matched$expr, GENE_CLASS, FIGURES_DIR_01)
  violin_p <- plot_expression_violin(matched$expr, GENE_CLASS, FIGURES_DIR_01)
  heatmap_p <- plot_expression_heatmap(matched$expr, GENE_CLASS, FIGURES_DIR_01)

  # --- 6. Expression summary table ---
  expr_summary <- build_expression_summary_table(
    matched$expr, GENE_CLASS, found_genes, missing_genes
  )
  save_table_pub(expr_summary, file.path(TABLES_DIR_01, "Table01_expression_gene_summary"))

  # --- 7. Survival analysis ---
  cli::cli_h1("Survival Analysis")

  km_results <- list()
  skipped_genes <- c()

  for (gene in found_genes) {
    cli::cli_h3("Analyzing {.val {gene}}")

    if (matched$n_matched < MIN_PATIENTS_KM) {
      cli::cli_alert_warning("Skipping: only {.val {matched$n_matched}} patients (< {.val {MIN_PATIENTS_KM}})")
      skipped_genes <- c(skipped_genes, gene)
      next
    }

    expr_vec <- matched$expr[gene, ]
    km_res <- run_km_for_gene(gene, expr_vec, matched$clinical)

    # Check minimum events
    if (km_res$n_events < 2 * MIN_EVENTS_PER_GROUP) {
      cli::cli_alert_warning("Skipping: only {.val {km_res$n_events}} events (< {.val {2 * MIN_EVENTS_PER_GROUP}} required)")
      skipped_genes <- c(skipped_genes, gene)
      next
    }

    km_results[[gene]] <- km_res

    # Individual KM plot
    plot_km_gene(km_res, GENE_CLASS, FIGURES_DIR_01)

    cli::cli_alert_info(
      "HR={.val {round(km_res$cox_HR, 2)}} [{.val {round(km_res$cox_CI_lower, 2)}}–{.val {round(km_res$cox_CI_upper, 2)}}], ",
      "log-rank P={.val {format_pval(km_res$logrank_p)}}, ",
      "n={.val {km_res$n_patients}}, events={.val {km_res$n_events}}"
    )
  }

  # --- 8. Survival summary table ---
  if (length(km_results) > 0) {
    surv_summary <- build_survival_summary_table(km_results)
    save_table_pub(surv_summary, file.path(TABLES_DIR_01, "Table02_survival_summary"))

    # Combined KM figure
    plot_km_combined(km_results, GENE_CLASS, FIGURES_DIR_01)

    # Forest plot
    plot_survival_forest(surv_summary, FIGURES_DIR_01)

    # --- Phase 2: robustness analyses (continuous / multivariable / cutpoint) ---
    adv_summary <- run_advanced_survival(km_results, matched$clinical, surv_summary)
    save_table_pub(adv_summary, file.path(TABLES_DIR_01, "Table04_survival_robustness"))
    plot_continuous_cox_forest(adv_summary, FIGURES_DIR_01)
    write_robustness_notes(adv_summary, OUTPUT_DIR_01)
  } else {
    cli::cli_alert_warning("No genes passed survival analysis thresholds")
    surv_summary <- data.frame()
    adv_summary <- data.frame()
  }

  # --- 9. Sample QC table ---
  sample_qc <- build_sample_qc_table(
    n_expression_aliquots = ncol(se),
    n_unique_patients = ncol(expr_tumor_patient),
    n_tumor_samples = sum(sample_info$sample_type == "Primary Tumor"),
    n_normal_samples = sum(sample_info$sample_type == "Solid Tissue Normal"),
    n_matched_clinical = matched$n_matched,
    found_genes = found_genes,
    missing_genes = missing_genes
  )
  save_table_pub(sample_qc, file.path(TABLES_DIR_01, "Table03_sample_qc_summary"))

  # --- 10. Module summary ---
  write_module_summary(
    expr_summary, surv_summary, sample_qc,
    found_genes, missing_genes,
    matched$n_matched,
    sum(sample_info$sample_type == "Solid Tissue Normal"),
    OUTPUT_DIR_01
  )

  # --- 11. Session info ---
  export_session_info(OUTPUT_DIR_01)

  # --- Done ---
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cli::cli_h1("Module 01 Complete")
  cli::cli_alert_success("Finished in {.val {round(elapsed, 1)}} minutes")
  cli::cli_alert_info("Outputs in {.path {OUTPUT_DIR_01}}")
  cli::cli_bullets(c(
    "v" = sprintf("Expression: %d genes analyzed (%d missing)", length(found_genes), length(missing_genes)),
    "v" = sprintf("Survival: %d genes passed QC (%d KM analyses)", length(km_results), length(km_results)),
    "v" = sprintf("Figures: %s", FIGURES_DIR_01),
    "v" = sprintf("Tables: %s", TABLES_DIR_01)
  ))

  invisible(list(
    se = se,
    expr_target = expr_target,
    clinical = clinical,
    matched = matched,
    km_results = km_results,
    surv_summary = surv_summary,
    adv_summary = adv_summary,
    expr_summary = expr_summary,
    sample_qc = sample_qc,
    elapsed = elapsed
  ))
}

# ============================================================================
# Execute
# ============================================================================

results <- run_pipeline()

# Keep objects available for interactive inspection
message("\n[01_expression_survival] Pipeline complete. Results stored in 'results' list.")
