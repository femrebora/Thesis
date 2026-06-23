library(here)
###############################################################################
# Metabolomics Figures — A172-S vs A172-R
# Glioblastoma metabolic profiling
#
# This script regenerates publication-quality figures for:
#   A172-S (Sensitive, n=4) vs A172-R (Resistant, n=4)
#
# Dependencies: MetaboAnalystR 4.0.0, ggplot2, pheatmap, EnhancedVolcano
# R version: 4.4.2
###############################################################################

# ============================================================================
# 0. SETUP
# ============================================================================
library(MetaboAnalystR)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(tidyverse)
library(mixOmics)
library(ggrepel)

# Project paths
PROJ_DIR <- here::here("metabolomics")
DATA_DIR   <- file.path(PROJ_DIR, "data")
FIG_DIR    <- file.path(PROJ_DIR, "figures")
SCRIPT_DIR <- file.path(PROJ_DIR, "scripts")



# Color scheme (ColorBrewer Set1, colorblind-friendly)
COL_R <- "#E41A1C"  # A172-R (Resistant) — red
COL_S <- "#377EB8"  # A172-S (Sensitive) — blue
GROUP_COLORS <- c("A172-R" = COL_R, "A172-S" = COL_S)

# High-resolution output settings
DPI         <- 300
FIGURE_W    <- 8    # inches
FIGURE_H    <- 7    # inches

# Helper: save figure in multiple formats
save_fig <- function(p, basename, w = FIGURE_W, h = FIGURE_H) {
  ggsave(file.path(FIG_DIR, paste0(basename, ".pdf")), plot = p,
         width = w, height = h, device = "pdf")
  ggsave(file.path(FIG_DIR, paste0(basename, ".png")), plot = p,
         width = w, height = h, dpi = DPI)
  ggsave(file.path(FIG_DIR, paste0(basename, ".tiff")), plot = p,
         width = w, height = h, dpi = DPI, device = "tiff",
         compression = "lzw")
}

# ============================================================================
# 1. DATA LOADING & PREPROCESSING
# ============================================================================
cat("\n=== Phase 1: Data Loading ===\n")

# Initialize MetaboAnalystR for peak intensity table (statistical analysis)
mSet <- InitDataObjects("pktable", "stat", FALSE)

# Read the cleaned data
input_file <- file.path(DATA_DIR, "metabolomics_clean.csv")
mSet <- Read.TextData(mSet, input_file, "colu", "disc")

cat("Data dimensions:", dim(mSet$dataSet$orig), "\n")
cat("Groups found:", unique(mSet$dataSet$cls), "\n")

# --- Sanity checks ---
mSet <- SanityCheckData(mSet)

# --- Missing value imputation ---
# Replace zeros / min values with 1/5 of the minimum positive value
# (standard practice for LC-MS metabolomics)
mSet <- ReplaceMin(mSet)

# --- Filtering ---
# Remove features with >25% missing values, using IQR filter
mSet <- FilterVariable(mSet, "none", 25, "none", -1, "mean", 0)
cat("Features after filtering:", nrow(mSet$dataSet$proc), "\n")

# --- Normalization ---
# Strategy: try multiple normalization approaches and compare
# Final choice: MedianNorm + LogNorm + ParetoNorm (standard for metabolomics)

cat("\n=== Phase 2: Normalization ===\n")

# Apply standard normalization: Median + Log + Pareto
# This is the gold standard for metabolomics (reduces heteroscedasticity
# while preserving biological variation)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm",
                      ratio = FALSE, ratioNum = 20)

# Save normalized data
norm_data <- as.data.frame(t(mSet$dataSet$norm))
write.csv(norm_data, file.path(DATA_DIR, "data_normalized_final.csv"))
cat("Normalized data saved to data_normalized_final.csv\n")

# ============================================================================
# 2. STATISTICAL ANALYSIS
# ============================================================================
cat("\n=== Phase 3: Statistical Analysis ===\n")

# --- Fold Change Analysis ---
mSet <- FC.Anal(mSet, fc.thresh = 2.0, cmp.type = 0, paired = FALSE)

# Extract fold change data for custom plots (v4.0.0: list of named vectors)
fc_log2 <- mSet$analSet$fc$fc.log   # log2 fold change
fc_ratio <- mSet$analSet$fc$fc.all  # raw fold change ratio
metab_names_fc <- names(fc_log2)

fc_data <- data.frame(
  metabolite = metab_names_fc,
  log2FC = fc_log2,
  FC = fc_ratio,
  row.names = NULL,
  stringsAsFactors = FALSE
)

# --- T-tests ---
# Welch's t-test (unequal variance) and Student's t-test (equal variance)
mSet <- Ttests.Anal(mSet, nonpar = FALSE, threshp = 0.05, paired = FALSE,
                    equal.var = FALSE, pvalType = "fdr")
mSet <- Ttests.Anal(mSet, nonpar = FALSE, threshp = 0.05, paired = FALSE,
                    equal.var = TRUE, pvalType = "fdr")

# Extract t-test results (v4.0.0: list of named vectors)
tt_data <- data.frame(
  metabolite  = names(mSet$analSet$tt$t.score),
  t_score     = mSet$analSet$tt$t.score,
  p_value     = mSet$analSet$tt$p.value,
  fdr_p       = mSet$analSet$tt$fdr.p,
  neg_log10_p = mSet$analSet$tt$p.log,
  significant = mSet$analSet$tt$inx.imp,
  stringsAsFactors = FALSE
)

cat("T-test results: number of features with raw p < 0.05:",
    sum(tt_data$p_value < 0.05, na.rm = TRUE), "\n")
cat("T-test results: number of features with FDR < 0.05:",
    sum(tt_data$fdr_p < 0.05, na.rm = TRUE), "\n")

# --- Wilcoxon rank-sum test ---
mSet <- Ttests.Anal(mSet, nonpar = TRUE, threshp = 0.05, paired = FALSE,
                    equal.var = FALSE, pvalType = "fdr")

# --- Volcano plot (MetaboAnalystR built-in) ---
mSet <- Volcano.Anal(mSet, paired = FALSE, fcthresh = 2.0, cmpType = 0,
                     nonpar = FALSE, threshp = 0.05, equal.var = FALSE,
                     pval.type = "raw")

# ============================================================================
# 3. MULTIVARIATE ANALYSIS (base R + mixOmics)
# ============================================================================
cat("\n=== Phase 4: Multivariate Analysis ===\n")

# Get the normalized data from MetaboAnalystR
# v4.0.0: mSet$dataSet$norm is already samples (rows) × metabolites (columns)
# Sample order: rows 1-4 = Resistant (R1-R4), rows 5-8 = Sensitive (S1-S4)
norm_data <- mSet$dataSet$norm  # data.frame, samples × metabolites
norm_data <- as.data.frame(norm_data)  # ensure it's a data.frame
sample_order <- c("R1", "R2", "R3", "R4", "S1", "S2", "S3", "S4")
sample_groups <- c(rep("A172-R", 4), rep("A172-S", 4))
rownames(norm_data) <- sample_order

# --- PCA with prcomp ---
cat("Running PCA (prcomp)...\n")
pca_res <- prcomp(norm_data, center = TRUE, scale. = FALSE)
# Note: data already Pareto-scaled, so scale.=FALSE

pca_var <- round(100 * summary(pca_res)$importance[2, ], 1)

pca_scores <- as.data.frame(pca_res$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Sample <- rownames(pca_scores)
pca_scores$Group <- sample_groups
pca_scores$Label <- sample_order

# --- PLS-DA with mixOmics ---
cat("Running PLS-DA (mixOmics)...\n")

# Prepare data for mixOmics: rows = samples, columns = variables
X <- as.matrix(norm_data)  # already 8 samples × 104 metabolites
Y <- factor(sample_groups)

# Fit PLS-DA with 2 components
plsda_res <- plsda(X, Y, ncomp = 2)

# Extract scores
plsda_scores <- as.data.frame(plsda_res$variates$X)
colnames(plsda_scores) <- c("Comp1", "Comp2")
plsda_scores$Sample <- rownames(plsda_scores)
plsda_scores$Group <- sample_groups
plsda_scores$Label <- sample_order

# PLS-DA variance explained
plsda_var_x <- round(100 * plsda_res$prop_expl_var$X, 1)

# --- PLS-DA VIP scores ---
# VIP is calculated as weighted sum of squared loadings
plsda_loadings <- plsda_res$loadings$X
plsda_w <- plsda_loadings[, 1]  # weights for component 1
plsda_w2 <- plsda_loadings[, 2] # weights for component 2

# Calculate VIP for comp 1 + comp 2 combined
ssx <- apply(plsda_res$variates$X, 2, function(v) sum(v^2))
vip_scores <- sqrt(ncol(X) * (ssx[1] * plsda_w^2 + ssx[2] * plsda_w2^2) /
                   sum(ssx))
vip_scores <- sort(vip_scores, decreasing = TRUE)

plsda_vip <- data.frame(
  metabolite = names(vip_scores),
  VIP = as.numeric(vip_scores),
  stringsAsFactors = FALSE
)

# --- PLS-DA Cross-Validation ---
cat("Running PLS-DA cross-validation...\n")
perf_plsda <- perf(plsda_res, validation = "loo", progressBar = FALSE)
cat("PLS-DA performance:\n")
print(perf_plsda)

# --- Custom OPLS-DA-like analysis ---
# For 2-group comparison, PLS-DA with 1 predictive + 1 orthogonal component
# gives equivalent information to OPLS-DA. We compute loadings and
# use them to identify discriminative metabolites.

# Compute S-plot-like metrics: correlation between X and the PLS-DA scores
plsda_t1 <- plsda_res$variates$X[, 1]
splot_df <- data.frame(
  metabolite = colnames(X),
  covariance = as.numeric(cor(X, plsda_t1, use = "complete.obs")),
  row.names = NULL,
  stringsAsFactors = FALSE
)
splot_df$correlation <- abs(splot_df$covariance)
splot_df$VIP <- plsda_vip$VIP[match(splot_df$metabolite, plsda_vip$metabolite)]
splot_df$label <- ifelse(splot_df$VIP > 1.3 | abs(splot_df$covariance) > 0.6,
                         splot_df$metabolite, "")

cat("Top 5 metabolites by VIP:\n")
print(head(plsda_vip, 5))

# ============================================================================
# 4. PUBLICATION-QUALITY FIGURES
# ============================================================================
cat("\n=== Phase 5: Generating Publication-Quality Figures ===\n")

# Common theme for publication-quality figures
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
  )

# --------------------------------------------------------------------------
# Figure 1: PCA Score Plot
# --------------------------------------------------------------------------
cat("  Figure 1: PCA Score Plot\n")

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  stat_ellipse(level = 0.95, linewidth = 1.0) +
  geom_point(size = 4.5, alpha = 0.85) +
  geom_text(aes(label = Label), vjust = -1.2, size = 3.5, show.legend = FALSE, fontface = "bold") +
  scale_color_manual(values = GROUP_COLORS) +
  labs(
    title = "PCA Score Plot",
    subtitle = paste0("A172-S vs A172-R (", nrow(pca_scores), " samples, ",
                      nrow(mSet$dataSet$norm), " metabolites)"),
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)"),
    color = "Group"
  ) +
  theme_pub +
  coord_fixed(ratio = 1)
save_fig(p_pca, "pca_score2d")

# PCA Scree Plot
pca_all_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
scree_df <- data.frame(
  PC = factor(paste0("PC", 1:length(pca_all_var)),
              levels = paste0("PC", 1:length(pca_all_var))),
  Variance = pca_all_var
)

p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = COL_S, alpha = 0.8) +
  geom_line(aes(group = 1), linewidth = 0.8, color = "red3") +
  geom_point(size = 2, color = "red3") +
  labs(title = "PCA Scree Plot", x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_scree, "pca_scree", w = 10)

# --------------------------------------------------------------------------
# Figure 2: PLS-DA Score Plot
# --------------------------------------------------------------------------
cat("  Figure 2: PLS-DA Score Plot\n")

# Get PLS-DA variance explained from model
plsda_var_x <- round(100 * mSet$analSet$plsr$Xvar / sum(mSet$analSet$plsr$Xvar), 1)

p_plsda <- ggplot(plsda_scores, aes(x = Comp1, y = Comp2, color = Group)) +
  stat_ellipse(level = 0.95, linewidth = 1.0) +
  geom_point(size = 4.5, alpha = 0.85) +
  geom_text(aes(label = Label), vjust = -1.2, size = 3.5, show.legend = FALSE, fontface = "bold") +
  scale_color_manual(values = GROUP_COLORS) +
  labs(
    title = "PLS-DA Score Plot",
    subtitle = "A172-S vs A172-R",
    x = paste0("Component 1 (", plsda_var_x[1], "%)"),
    y = paste0("Component 2 (", plsda_var_x[2], "%)"),
    color = "Group"
  ) +
  theme_pub +
  coord_fixed(ratio = 1)
save_fig(p_plsda, "plsda_score2d")

# --------------------------------------------------------------------------
# Figure 3: PLS-DA VIP Scores (Top 20)
# --------------------------------------------------------------------------
cat("  Figure 3: PLS-DA VIP Plot\n")

vip_top20 <- head(plsda_vip, 20)
vip_top20$metabolite <- factor(vip_top20$metabolite,
                               levels = rev(vip_top20$metabolite))

p_vip <- ggplot(vip_top20, aes(x = VIP, y = metabolite)) +
  geom_bar(stat = "identity", fill = COL_S, alpha = 0.85) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "red3", linewidth = 0.8) +
  annotate("text", x = 1.05, y = 1, label = "VIP = 1",
           hjust = 0, size = 3.5, color = "red3") +
  labs(title = "PLS-DA Variable Importance in Projection (VIP)",
       subtitle = "Top 20 metabolites", x = "VIP Score", y = "") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9))
save_fig(p_vip, "plsda_vip", w = 10, h = 8)

# --------------------------------------------------------------------------
# Figure 4: Cross-validated PLS-DA Performance
# --------------------------------------------------------------------------
cat("  Figure 4: PLS-DA CV Performance\n")

# Plot cross-validation results
pdf(file.path(FIG_DIR, "plsda_cv_performance.pdf"), width = 10, height = 7)
plot(perf_plsda, overlay = "measure", legend = FALSE)
dev.off()

png(file.path(FIG_DIR, "plsda_cv_performance.png"), width = 10, height = 7,
    units = "in", res = DPI)
plot(perf_plsda, overlay = "measure", legend = FALSE)
dev.off()
cat("  PLS-DA CV performance saved\n")

# --------------------------------------------------------------------------
# Figure 5: PLS-DA S-Plot (Correlation-based)
# --------------------------------------------------------------------------
cat("  Figure 5: PLS-DA S-Plot\n")

p_splot <- ggplot(splot_df, aes(x = covariance, y = correlation)) +
  geom_point(aes(size = VIP, color = VIP), alpha = 0.7) +
  scale_color_gradient(low = "grey60", high = "red3") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15,
                  box.padding = 0.3) +
  labs(title = "PLS-DA S-Plot",
       subtitle = "Correlation with Component 1 vs Covariance (colored by VIP)",
       x = "Covariance with Comp 1", y = "|Correlation| with Comp 1") +
  theme_pub +
  theme(legend.position = "right")
save_fig(p_splot, "plsda_splot", w = 10)

# --------------------------------------------------------------------------
# Figure 6: Volcano Plot
# --------------------------------------------------------------------------
cat("  Figure 6: Volcano Plot\n")

# Use t-test results for volcano (Welch's t-test p-values)
volc_df <- tt_data
volc_df$log2FC <- fc_data$log2FC[match(volc_df$metabolite, fc_data$metabolite)]
volc_df$neg_log10_p <- volc_df$neg_log10_p  # already computed

# Define significance thresholds
volc_df$significance <- "Not Significant"
volc_df$significance[volc_df$p_value < 0.05 & abs(volc_df$log2FC) > 1] <- "p < 0.05 & |log2FC| > 1"
volc_df$significance[volc_df$p_value < 0.05] <- "p < 0.05"

# Label top metabolites by combined ranking
volc_df$rank_score <- abs(volc_df$log2FC) * volc_df$neg_log10_p
volc_df <- volc_df[order(-volc_df$rank_score), ]
volc_df$label <- ""
volc_df$label[1:15] <- volc_df$metabolite[1:15]

sig_colors <- c("Not Significant" = "grey70",
                "p < 0.05" = "orange",
                "p < 0.05 & |log2FC| > 1" = "red3")

p_volcano <- ggplot(volc_df, aes(x = log2FC, y = neg_log10_p,
                                  color = significance, label = label)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = sig_colors) +
  geom_text(size = 2.8, vjust = -0.8, show.legend = FALSE, check_overlap = TRUE,
            color = "black") +
  labs(title = "Volcano Plot: A172-R vs A172-S",
       subtitle = paste0("NOTE: No metabolites survive FDR correction (n=4/group). ",
                         "Raw p-values shown for exploratory use."),
       x = "log2(Fold Change) [A172-R / A172-S]",
       y = "-log10(p-value) [Welch's t-test]",
       color = "Significance",
       caption = "Dashed lines: |log2FC| = 1, p = 0.05 (unadjusted)") +
  theme_pub +
  theme(plot.caption = element_text(size = 9, color = "grey50", hjust = 0))
save_fig(p_volcano, "volcano", w = 10, h = 8)

# --------------------------------------------------------------------------
# Also create an EnhancedVolcano version
# --------------------------------------------------------------------------
cat("  Figure 6b: Enhanced Volcano Plot\n")

enhanced_volc_df <- data.frame(
  gene = volc_df$metabolite,
  log2FC = volc_df$log2FC,
  pvalue = volc_df$p_value,
  row.names = volc_df$metabolite
)

pdf(file.path(FIG_DIR, "volcano_enhanced.pdf"), width = 10, height = 9)
EnhancedVolcano(enhanced_volc_df,
                lab = rownames(enhanced_volc_df),
                x = "log2FC", y = "pvalue",
                title = "A172-R vs A172-S",
                subtitle = "No FDR correction applied (n=4/group)",
                pCutoff = 0.05, FCcutoff = 1.0,
                pointSize = 3.0, labSize = 3.5,
                col = c("grey30", "grey30", "green3", "red2"),
                colAlpha = 0.6,
                legendPosition = "right",
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                max.overlaps = 20,
                caption = paste0("Total: ", nrow(enhanced_volc_df), " metabolites"))
dev.off()
cat("  Enhanced volcano saved to volcano_enhanced.pdf\n")

# --------------------------------------------------------------------------
# Figure 7: Heatmap of Top Differential Metabolites
# --------------------------------------------------------------------------
cat("  Figure 7: Heatmap\n")

# Create combined ranking: VIP (from PLS-DA) + abs(log2FC) + -log10(p)
combined_rank <- plsda_vip
combined_rank$log2FC <- fc_data$log2FC[match(combined_rank$metabolite, fc_data$metabolite)]
combined_rank$p_value <- tt_data$p_value[match(combined_rank$metabolite, tt_data$metabolite)]
combined_rank$neg_log10_p <- -log10(combined_rank$p_value)

# Normalize each metric to 0-1 and combine
combined_rank$vip_norm <- combined_rank$VIP / max(combined_rank$VIP, na.rm = TRUE)
combined_rank$fc_norm <- abs(combined_rank$log2FC) / max(abs(combined_rank$log2FC), na.rm = TRUE)
combined_rank$p_norm <- combined_rank$neg_log10_p / max(combined_rank$neg_log10_p, na.rm = TRUE)
combined_rank$combined_score <- with(combined_rank,
                                     vip_norm + fc_norm + p_norm)
combined_rank <- combined_rank[order(-combined_rank$combined_score), ]

# Select top N metabolites for heatmap
N_HEATMAP <- min(40, nrow(combined_rank))
top_metabolites <- combined_rank$metabolite[1:N_HEATMAP]

# Get normalized expression data for these metabolites
norm_mat <- as.data.frame(t(mSet$dataSet$norm))
# norm_data is samples × metabolites; subset columns and transpose to metabolites × samples
heatmap_data <- as.data.frame(t(norm_data[, top_metabolites, drop = FALSE]))
colnames(heatmap_data) <- sample_order

# Z-score rows
heatmap_data_z <- t(scale(t(heatmap_data)))

# Sample annotation
sample_groups_df <- data.frame(
  Group = sample_groups,
  row.names = sample_order
)
ann_colors <- list(Group = c("A172-R" = COL_R, "A172-S" = COL_S))

# Clean up very long metabolite names for display
rownames(heatmap_data_z) <- substr(rownames(heatmap_data_z), 1, 45)

# Generate heatmap
pdf(file.path(FIG_DIR, "heatmap.pdf"), width = 10, height = 12)
pheatmap(heatmap_data_z,
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         annotation_col = sample_groups_df,
         annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         main = paste0("Top ", N_HEATMAP, " Differential Metabolites\nA172-S vs A172-R"),
         fontsize_row = 8,
         fontsize = 10,
         border_color = NA,
         show_colnames = TRUE)
dev.off()

# Also PNG version
png(file.path(FIG_DIR, "heatmap.png"), width = 10, height = 12, units = "in", res = DPI)
pheatmap(heatmap_data_z,
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         annotation_col = sample_groups_df,
         annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         main = paste0("Top ", N_HEATMAP, " Differential Metabolites\nA172-S vs A172-R"),
         fontsize_row = 8,
         fontsize = 10,
         border_color = NA,
         show_colnames = TRUE)
dev.off()

# Also TIFF version
tiff(file.path(FIG_DIR, "heatmap.tiff"), width = 10, height = 12, units = "in",
     res = DPI, compression = "lzw")
pheatmap(heatmap_data_z,
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         annotation_col = sample_groups_df,
         annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         main = paste0("Top ", N_HEATMAP, " Differential Metabolites\nA172-S vs A172-R"),
         fontsize_row = 8,
         fontsize = 10,
         border_color = NA,
         show_colnames = TRUE)
dev.off()
cat("  Heatmap saved\n")

# --------------------------------------------------------------------------
# Figure 8: Boxplots of Top Metabolites
# --------------------------------------------------------------------------
cat("  Figure 8: Boxplots\n")

# Get original (normalized) data for top metabolites
N_BOX <- 9
top_for_box <- combined_rank$metabolite[1:N_BOX]

# Build tidy data frame for boxplots
box_list <- list()
for (metab in top_for_box) {
  vals <- as.numeric(norm_data[, metab])
  samples <- sample_order
  groups <- sample_groups
  box_list[[metab]] <- data.frame(
    Metabolite = metab,
    Sample = samples,
    Group = groups,
    Intensity = vals,
    stringsAsFactors = FALSE
  )
}
box_df <- do.call(rbind, box_list)

# Shorten long metabolite names
box_df$Metabolite_short <- substr(box_df$Metabolite, 1, 40)

# Add p-values
pv_df <- data.frame(
  Metabolite = top_for_box,
  p_value = combined_rank$p_value[match(top_for_box, combined_rank$metabolite)]
)
pv_df$Metabolite_short <- substr(pv_df$Metabolite, 1, 40)
pv_df$p_label <- paste0("p = ", formatC(pv_df$p_value, format = "e", digits = 2))

box_df$Metabolite_short <- factor(box_df$Metabolite_short,
                                   levels = unique(pv_df$Metabolite_short))

p_box <- ggplot(box_df, aes(x = Group, y = Intensity, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.9, color = "black") +
  geom_text(data = pv_df, aes(x = 1.5, y = Inf, label = p_label),
            vjust = 1.5, size = 3.0, inherit.aes = FALSE, color = "grey40") +
  scale_fill_manual(values = GROUP_COLORS) +
  facet_wrap(~ Metabolite_short, scales = "free_y", ncol = 3) +
  labs(title = "Top Differential Metabolites: A172-R vs A172-S",
       subtitle = paste0("Selected by combined VIP + fold-change + p-value ranking. ",
                         "Raw Welch's t-test p-values shown (no FDR correction)."),
       x = "", y = "Normalized Intensity (Log, Pareto-scaled)") +
  theme_pub +
  theme(strip.text = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 10))
save_fig(p_box, "boxplots_top9", w = 14, h = 12)

# --------------------------------------------------------------------------
# Figure 9: PLS-DA VIP Top 20
# --------------------------------------------------------------------------
cat("  Figure 9: PLS-DA VIP\n")

plsda_vip_top20 <- head(plsda_vip, 20)
plsda_vip_top20$metabolite <- factor(plsda_vip_top20$metabolite,
                                      levels = rev(plsda_vip_top20$metabolite))

p_vip2 <- ggplot(plsda_vip_top20, aes(x = VIP, y = metabolite)) +
  geom_bar(stat = "identity",
           fill = ifelse(plsda_vip_top20$VIP > 1, COL_R, "steelblue"),
           alpha = 0.85) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "red3", linewidth = 0.8) +
  labs(title = "PLS-DA Variable Importance in Projection (VIP)",
       subtitle = "Top 20 metabolites", x = "VIP Score", y = "") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9))
save_fig(p_vip2, "plsda_vip_top20", w = 10, h = 8)

# ============================================================================
# 5. SESSION INFO & SAVE
# ============================================================================
cat("\n=== Phase 6: Saving Results ===\n")

# Save combined ranking table
write.csv(combined_rank, file.path(DATA_DIR, "combined_metabolite_ranking.csv"),
          row.names = FALSE)
cat("Combined ranking saved\n")

# Save PLS-DA VIP
write.csv(plsda_vip, file.path(DATA_DIR, "plsda_vip_scores.csv"), row.names = FALSE)
cat("PLS-DA VIP scores saved\n")

# Save t-test results with fold change
tt_export <- tt_data
tt_export$log2FC <- fc_data$log2FC[match(tt_export$metabolite, fc_data$metabolite)]
write.csv(tt_export, file.path(DATA_DIR, "ttest_with_foldchange.csv"), row.names = FALSE)
cat("T-test + fold change results saved\n")

# Session info
cat("\n=== Session Info ===\n")
sessionInfo()

cat("\n=== FIGURE GENERATION COMPLETE ===\n")
cat("All figures saved to:", FIG_DIR, "\n")
cat("All data saved to:", DATA_DIR, "\n")
