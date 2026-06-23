library(here)
###############################################################################
# Multi-Omics Integration Figures: CRISPR + Metabolomics
# A172 GBM: Sensitive vs Resistant
#
# Integrates:
#   - CRISPR-Seq: Sabatini metabolic gene library (15 sig genes, Spost_vs_Rpost)
#   - Metabolomics: LC-MS (104 metabolites, A172-S vs A172-R)
###############################################################################

library(ggplot2)
library(pheatmap)
library(tidyverse)

PROJ_DIR <- here::here("metabolomics")
FIG_DIR  <- file.path(PROJ_DIR, "figures")
DATA_DIR <- file.path(PROJ_DIR, "data")
DPI      <- 300



save_fig <- function(p, basename, w = 10, h = 8) {
  ggsave(file.path(FIG_DIR, paste0(basename, ".pdf")), plot = p, width = w, height = h, device = "pdf")
  ggsave(file.path(FIG_DIR, paste0(basename, ".png")), plot = p, width = w, height = h, dpi = DPI)
  ggsave(file.path(FIG_DIR, paste0(basename, ".tiff")), plot = p, width = w, height = h, dpi = DPI, device = "tiff", compression = "lzw")
}

theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "top",
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"))

COL_R <- "#E41A1C"  # A172-R
COL_S <- "#377EB8"  # A172-S

# ============================================================================
# 1. LOAD INTEGRATION DATA
# ============================================================================
cat("Loading integration data...\n")

# Read pathway integration table
int_df <- read.csv(file.path(DATA_DIR, "integrated_pathway_analysis.csv"))
cat("Pathways:", nrow(int_df), "\n")

# Read gene-metabolite connections
conn_df <- read.csv(file.path(DATA_DIR, "gene_metabolite_connections.csv"))
cat("Connections:", nrow(conn_df), "\n")

# Read individual datasets
tt_df <- read.csv(file.path(DATA_DIR, "ttest_with_foldchange.csv"))

# CRISPR genes with manual data
crispr_genes <- data.frame(
  Gene = c("PLA2G4E","SLC1A1","RRM2","ARSD","GRIN1","PFKFB1","SLC25A15",
           "SLC6A3","NAT2","TECR","TKTL1","APOC1","CYB5A","SLC25A20","ALOX15"),
  LFC = c(-12.101, -11.252, -10.764, -10.209, -10.181, -9.620, -9.458,
          -9.151, 1.615, -0.814, -0.292, -0.292, 0.284, 1.230, 2.373),
  pval = c(0.0026, 0.0083, 0.0145, 0.0208, 0.0265, 0.0324, 0.0383,
           0.0443, 0.0061, 0.0408, 0.0352, 0.0352, 0.0265, 0.0145, 0.0083),
  FDR = c(0.428, 0.683, 0.792, 0.851, 0.870, 0.886, 0.898,
          0.907, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636, 0.636),
  Direction = c(rep("Depleted", 8), rep("Enriched", 7)),
  stringsAsFactors = FALSE
)
crispr_genes$absLFC <- abs(crispr_genes$LFC)
crispr_genes$neg_log10_p <- -log10(crispr_genes$pval)

# ============================================================================
# 2. FIGURE 1: DUAL VOLCANO — CRISPR + Metabolomics
# ============================================================================
cat("Figure 1: Dual Volcano Plot\n")

# Prepare CRISPR volcano data
crispr_genes$omic <- "CRISPR"
crispr_genes$label <- crispr_genes$Gene

# Prepare metabolomics volcano data
metab_volc <- tt_df
metab_volc$omic <- "Metabolomics"
metab_volc$label <- metab_volc$metabolite
# Truncate long names
metab_volc$label <- ifelse(nchar(metab_volc$label) > 25,
                           paste0(substr(metab_volc$label, 1, 22), "..."),
                           metab_volc$label)

# Mark significance
crispr_genes$sig <- ifelse(crispr_genes$FDR < 0.9, "FDR < 0.9", "Not sig")
metab_volc$sig <- ifelse(metab_volc$fdr_p < 0.05, "FDR < 0.05",
                         ifelse(metab_volc$p_value < 0.05, "p < 0.05", "Not sig"))

# Select top labels (not all — too crowded)
crispr_top <- crispr_genes[1:10, ]  # top 10 by p-value
metab_top <- metab_volc[metab_volc$p_value < 0.05, ]
metab_top <- metab_top[order(metab_top$p_value), ]
metab_top <- head(metab_top, 10)

p_dual <- ggplot() +
  # CRISPR layer
  geom_point(data = crispr_genes, aes(x = LFC, y = neg_log10_p, color = "CRISPR"), size = 3, alpha = 0.8) +
  geom_text(data = crispr_top, aes(x = LFC, y = neg_log10_p, label = Gene),
            vjust = -0.8, size = 3, color = COL_R, fontface = "bold") +
  # Metabolomics layer
  geom_point(data = metab_volc, aes(x = log2FC, y = neg_log10_p, color = "Metabolomics"), size = 2.5, alpha = 0.7) +
  geom_text(data = metab_top, aes(x = log2FC, y = neg_log10_p, label = label),
            vjust = -0.8, size = 2.5, color = COL_S) +
  # Thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", alpha = 0.5) +
  # Styling
  scale_color_manual(values = c("CRISPR" = COL_R, "Metabolomics" = COL_S), name = "Omic Layer") +
  labs(title = "Multi-Omics Integration: CRISPR + Metabolomics",
       subtitle = "A172 Sensitive vs Resistant (Post-treatment)",
       x = "Log Fold Change (CRISPR: MAGeCK LFC  |  Metabolomics: log2FC)",
       y = "-log10(p-value)") +
  theme_pub +
  facet_wrap(~ omic, scales = "free_x", ncol = 2) +
  theme(strip.text = element_text(size = 12, face = "bold"))
save_fig(p_dual, "integration_dual_volcano", w = 16, h = 8)

# ============================================================================
# 3. FIGURE 2: PATHWAY OVERLAP HEATMAP
# ============================================================================
cat("Figure 2: Pathway Overlap Heatmap\n")

# Prepare pathway overlap matrix
# Rows: pathways with dual evidence
# Ensure correct types
int_df$Dual_Evidence <- as.logical(int_df$Dual_Evidence)
int_df$Combined_Score <- as.numeric(int_df$Combined_Score)
dual_pathways <- int_df[int_df$Dual_Evidence | (int_df$Combined_Score >= 2), ]

if (nrow(dual_pathways) >= 3) {
  # Build matrix: rows=pathways, cols=evidence type
  pw_mat <- as.matrix(dual_pathways[, c("CRISPR_Genes", "Metabolites")])
  rownames(pw_mat) <- dual_pathways$Pathway
  colnames(pw_mat) <- c("CRISPR Genes", "Metabolites")

  pdf(file.path(FIG_DIR, "integration_pathway_heatmap.pdf"), width = 8, height = max(6, nrow(pw_mat) * 0.4))
  pheatmap(pw_mat,
           color = colorRampPalette(c("white", "#FFF7BC", "#E41A1C"))(50),
           cluster_rows = TRUE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.0f",
           main = "Multi-Omics Pathway Evidence\n(CRISPR + Metabolomics)",
           fontsize = 10, fontsize_number = 10)
  dev.off()

  png(file.path(FIG_DIR, "integration_pathway_heatmap.png"), width = 8,
      height = max(6, nrow(pw_mat) * 0.4), units = "in", res = DPI)
  pheatmap(pw_mat,
           color = colorRampPalette(c("white", "#FFF7BC", "#E41A1C"))(50),
           cluster_rows = TRUE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.0f",
           main = "Multi-Omics Pathway Evidence\n(CRISPR + Metabolomics)",
           fontsize = 10, fontsize_number = 10)
  dev.off()

  tiff(file.path(FIG_DIR, "integration_pathway_heatmap.tiff"), width = 8,
       height = max(6, nrow(pw_mat) * 0.4), units = "in", res = DPI, compression = "lzw")
  pheatmap(pw_mat,
           color = colorRampPalette(c("white", "#FFF7BC", "#E41A1C"))(50),
           cluster_rows = TRUE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.0f",
           main = "Multi-Omics Pathway Evidence\n(CRISPR + Metabolomics)",
           fontsize = 10, fontsize_number = 10)
  dev.off()
  cat("  Pathway heatmap saved.\n")
}

# ============================================================================
# 4. FIGURE 3: GENE-METABOLITE NETWORK
# ============================================================================
cat("Figure 3: Gene-Metabolite Connection Plot\n")

if (nrow(conn_df) > 0) {
  # Fix column type
  conn_df$Gene_Direction <- as.character(conn_df$Gene_Direction)
  conn_df$Gene_Direction[is.na(conn_df$Gene_Direction)] <- "Depleted"

  p_conn <- ggplot(conn_df, aes(x = Gene, y = Metabolite, size = -log10(Metab_pvalue), color = Gene_Direction)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c("Depleted" = COL_R, "Enriched" = "#4DAF4A"), name = "CRISPR Direction") +
    scale_size_continuous(range = c(3, 8), name = "-log10(p)") +
    labs(title = "Gene-Metabolite Connections",
         subtitle = paste0("Shared metabolic pathways (", nrow(conn_df), " connections)"),
         x = "CRISPR Gene", y = "Metabolite") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9))
  save_fig(p_conn, "integration_gene_metab_network", w = 10, h = 6)
} else {
  cat("  No gene-metabolite connections to plot.\n")
}

# ============================================================================
# 5. FIGURE 4: CROSS-OMICS RANK PLOT
# ============================================================================
cat("Figure 4: Cross-Omics Rank Plot\n")

# Combine top genes and metabolites into a unified ranking
crispr_rank <- data.frame(
  Name = crispr_genes$Gene,
  Omic = "CRISPR",
  Rank = 1:15,
  LFC = crispr_genes$LFC,
  pval = crispr_genes$pval,
  stringsAsFactors = FALSE
)

metab_rank <- tt_df[order(tt_df$p_value), ]
metab_rank <- data.frame(
  Name = substr(metab_rank$metabolite[1:15], 1, 30),
  Omic = "Metabolomics",
  Rank = 1:15,
  LFC = metab_rank$log2FC[1:15],
  pval = metab_rank$p_value[1:15],
  stringsAsFactors = FALSE
)

combined_rank <- rbind(crispr_rank, metab_rank)
combined_rank$Omic <- factor(combined_rank$Omic, levels = c("CRISPR", "Metabolomics"))

p_rank <- ggplot(combined_rank, aes(x = Rank, y = abs(LFC), fill = Omic)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("CRISPR" = COL_R, "Metabolomics" = COL_S)) +
  labs(title = "Top Hits Comparison: CRISPR vs Metabolomics",
       subtitle = "Ranked by p-value, Top 15 in each omic layer",
       x = "Rank", y = "|Fold Change|") +
  theme_pub
save_fig(p_rank, "integration_rank_comparison", w = 12, h = 6)

# ============================================================================
# 6. FIGURE 5: PATHWAY COVERAGE DOT PLOT
# ============================================================================
cat("Figure 5: Pathway Coverage Dot Plot\n")

# Top pathways by combined score
top_pws <- head(int_df[order(-int_df$Combined_Score), ], 20)
top_pws$Pathway <- factor(top_pws$Pathway, levels = rev(top_pws$Pathway))

p_dot <- ggplot(top_pws, aes(x = CRISPR_Genes, y = Pathway, size = Metabolites, color = Combined_Score)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "steelblue", high = COL_R, name = "Combined Score") +
  scale_size_continuous(range = c(2, 8), name = "Metabolites") +
  geom_text(aes(label = paste0(CRISPR_Genes, "/", Metabolites)), hjust = -0.5, size = 3) +
  labs(title = "Multi-Omics Pathway Coverage",
       subtitle = paste0("CRISPR genes + Metabolites per pathway (", nrow(int_df), " pathways total)"),
       x = "Number of CRISPR Genes", y = "") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9))
save_fig(p_dot, "integration_pathway_dotplot", w = 12, h = 8)

cat("\n=== INTEGRATION FIGURES COMPLETE ===\n")
cat("All figures saved to:", FIG_DIR, "\n")

sessionInfo()
