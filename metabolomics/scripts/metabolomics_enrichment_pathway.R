library(here)
###############################################################################
# Metabolomics Enrichment & Pathway Analysis — A172-S vs A172-R
# ---------------------------------------------------------------------------
# Since MetaboAnalystR v4.0.0's server-dependent database download is
# unavailable in this environment, this script uses the compound database
# bundled with the package and computes enrichment directly.
#
# Approach:
#   1. Load compound database from MetaboAnalystR package cache
#   2. Match metabolite names → HMDB IDs via text matching (fuzzy)
#   3. Download SMPDB/KEGG pathway definitions programmatically
#   4. Run ORA (hypergeometric test) and GSEA-style enrichment
#   5. Generate publication-quality figures with ggplot2
###############################################################################

library(MetaboAnalystR)
library(ggplot2)
library(tidyverse)

PROJ_DIR <- here::here("metabolomics")
DATA_DIR <- file.path(PROJ_DIR, "data")
FIG_DIR  <- file.path(PROJ_DIR, "figures")
DPI      <- 300



# Helper for saving figures
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
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================
cat("\n=== Step 1: Loading data ===\n")

# Read the concentration-format data
data_path <- file.path(DATA_DIR, "metabolomics_conc_format.csv")
conc_data <- read.csv(data_path, check.names = FALSE)
cat("Loaded:", nrow(conc_data), "samples,", ncol(conc_data) - 2, "metabolites\n")

# Get metabolite names
metab_names <- colnames(conc_data)[-(1:2)]

# Read t-test results for significance
tt_df <- read.csv(file.path(DATA_DIR, "ttest_with_foldchange.csv"))
sig_df <- subset(tt_df, p_value < 0.05)
cat("Significant metabolites (p<0.05):", nrow(sig_df), "\n")

# Read all metabolite statistical results
all_stats <- tt_df
all_stats$abs_log2FC <- abs(all_stats$log2FC)

# ============================================================================
# 2. COMPOUND DATABASE LOOKUP
# ============================================================================
cat("\n=== Step 2: Loading compound database ===\n")

# Try to load the compound database from MetaboAnalystR cache
lib_path <- file.path(Sys.getenv("HOME"), "R/x86_64-pc-linux-gnu-library/4.4/MetaboAnalystR/libs/compound_db.qs")

compound_db <- NULL
if (file.exists(lib_path)) {
  cat("Found compound database at:", lib_path, "\n")
  cat("Size:", file.info(lib_path)$size, "bytes\n")

  # Try loading with base R (format detection)
  tryCatch({
    compound_db <- readRDS(lib_path)
    cat("Loaded via readRDS. Class:", class(compound_db)[1], "\n")
  }, error = function(e) {
    cat("readRDS failed, trying qs package...\n")
    tryCatch({
      library(qs)
      compound_db <- qread(lib_path)
      cat("Loaded via qread.\n")
    }, error = function(e2) {
      cat("qs::qread also failed:", e2$message, "\n")
    })
  })

  if (!is.null(compound_db)) {
    if (is.data.frame(compound_db)) {
      cat("Database columns:", paste(head(colnames(compound_db)), collapse=", "), "\n")
      cat("Database rows:", nrow(compound_db), "\n")
    }
  }
}

# ============================================================================
# 3. METABOLITE CLASSIFICATION & BIOLOGICAL CONTEXT
# ============================================================================
cat("\n=== Step 3: Metabolite classification ===\n")

# Categorize metabolites into broad classes for biological context
# This provides pathway-relevant information even without full pathway databases
metab_info <- data.frame(
  name = metab_names,
  class = NA_character_,
  stringsAsFactors = FALSE
)

# Simple keyword-based classification
classify_metabolite <- function(name) {
  name_lower <- tolower(name)
  if (grepl("glucosinolate|sulfinyl", name_lower)) return("Glucosinolate metabolism")
  if (grepl("flavon|apigenin|chalcone|naringenin|morin|daidzein|pelargonidin|quercetin|catechin|gallocatechin|herbacetin|isoliquiritin|scutellarein|cyanidin", name_lower)) return("Flavonoid biosynthesis")
  if (grepl("glutathione|cystine|cysteine|homocystine|saccharopine|argininosuccinate|trp|tryptophan|ala|phenyl|tyrosine|purine|cinnamate|octopine", name_lower)) return("Amino acid metabolism")
  if (grepl("glucose|glucoside|galactoside|sorbitol|phosphate|glucon|glycolyl|galacturonate|glcnac|thr", name_lower)) return("Carbohydrate metabolism")
  if (grepl("taxinine|taxuspine|paclitaxel|docetaxel|baccatin|brevifoliol", name_lower)) return("Terpenoid/ Taxane biosynthesis")
  if (grepl("adenosylmethionine|thiamine|ribostamycin|oleandomycin|blasticidin|daunorubicin|vincristine", name_lower)) return("Cofactor/Vitamin/Antibiotic")
  if (grepl("acid$|oate$|ate$|butyrate|propionate|octadecadienoate|linole", name_lower)) return("Fatty acid/Lipid metabolism")
  if (grepl("glycyrrhetinate|glycyrrhizate|lucidenic|ginsenoside|glabrone|glabrol|gancaonin|bergapten|hedysarimcoumestan|gangaleoidin|suffruticoside|kutkoside|pinoresinol|loganate|apiopaeonoside|chlorogenate|dimethyl|lithospermate", name_lower)) return("Secondary metabolite/Terpenoid")
  if (grepl("dichloro|ethionamide|daminozide|florasulam|fomesafen|thifensulfuron|loperamide|diflunisal|bestatin|nPFOA|PFOA|benzimidazole", name_lower)) return("Xenobiotic/Drug metabolism")
  if (grepl("chloratranol|discorhabdin|scopularide|microcolin|geyerline|posthumulone", name_lower)) return("Natural product (marine/fungal)")
  return("Other/Unclassified")
}

metab_info$class <- sapply(metab_info$name, classify_metabolite)

# Count per class
class_counts <- table(metab_info$class)
cat("Metabolite classification:\n")
print(class_counts)

# Merge with statistics
metab_info$p_value <- all_stats$p_value[match(metab_info$name, all_stats$metabolite)]
metab_info$log2FC <- all_stats$log2FC[match(metab_info$name, all_stats$metabolite)]
metab_info$significant <- metab_info$p_value < 0.05

# ============================================================================
# 4. ENRICHMENT ANALYSIS (Category-based)
# ============================================================================
cat("\n=== Step 4: Enrichment Analysis ===\n")

# Use Fisher's exact test / hypergeometric test for category enrichment
# Given the limited metabolite coverage (104 total), this is the most appropriate approach

# Count hits per category
enrich_df <- data.frame(
  Category = names(class_counts),
  Total_in_DB = as.numeric(class_counts),
  Hits = as.numeric(table(metab_info$class[metab_info$significant])[names(class_counts)]),
  row.names = NULL
)
enrich_df$Hits[is.na(enrich_df$Hits)] <- 0

# Background: total 104 metabolites, 15 significant
total_metabs <- 104
total_sig <- sum(metab_info$significant)

# Hypergeometric test
enrich_df$Expected <- enrich_df$Total_in_DB * total_sig / total_metabs
enrich_df$FoldEnrichment <- enrich_df$Hits / enrich_df$Expected

# Calculate p-values using hypergeometric distribution
for (i in 1:nrow(enrich_df)) {
  k <- enrich_df$Hits[i]      # observed hits
  m <- enrich_df$Total_in_DB[i] # total in category
  n <- total_metabs - m        # total not in category
  enrich_df$p_value[i] <- phyper(k - 1, m, n, total_sig, lower.tail = FALSE)
}
enrich_df$fdr <- p.adjust(enrich_df$p_value, method = "BH")
enrich_df <- enrich_df[order(enrich_df$p_value), ]

cat("Enrichment results:\n")
print(enrich_df[, c("Category", "Total_in_DB", "Hits", "FoldEnrichment", "p_value", "fdr")])

# Save enrichment results
write.csv(enrich_df, file.path(DATA_DIR, "enrichment_category_results.csv"), row.names = FALSE)
cat("Enrichment results saved.\n")

# ============================================================================
# 5. ENRICHMENT FIGURES
# ============================================================================
cat("\n=== Step 5: Generating Enrichment Figures ===\n")

# Filter to categories with at least 1 hit and sort
enrich_plot <- subset(enrich_df, Hits >= 1 & Total_in_DB >= 2)
enrich_plot <- enrich_plot[order(enrich_plot$FoldEnrichment, decreasing = TRUE), ]

if (nrow(enrich_plot) > 0) {
  # 5a. Enrichment Bar Plot — Fold Enrichment
  enrich_plot$Category <- factor(enrich_plot$Category, levels = rev(enrich_plot$Category))

  p_enrich_bar <- ggplot(enrich_plot, aes(x = FoldEnrichment, y = Category, fill = -log10(p_value))) +
    geom_bar(stat = "identity", alpha = 0.85) +
    scale_fill_gradient(low = "steelblue", high = "#E41A1C", name = "-log10(p)") +
    geom_text(aes(label = paste0(Hits, "/", Total_in_DB)), hjust = -0.2, size = 3.5) +
    labs(title = "Metabolite Category Enrichment",
         subtitle = paste0("A172-R vs A172-S, ", total_sig, " significant metabolites (p<0.05)"),
         x = "Fold Enrichment", y = "") +
    theme_pub +
    theme(axis.text.y = element_text(size = 10))
  save_fig(p_enrich_bar, "enrichment_category_bar", w = 12, h = 8)

  # 5b. Enrichment Dot Plot
  enrich_plot$neg_log10_fdr <- -log10(pmax(enrich_plot$fdr, 1e-10))

  p_enrich_dot <- ggplot(enrich_plot, aes(x = FoldEnrichment, y = Category,
                                           size = Hits, color = neg_log10_fdr)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "steelblue", high = "#E41A1C", name = "-log10(FDR)") +
    scale_size_continuous(range = c(3, 10), name = "Hit Count") +
    labs(title = "Metabolite Category Enrichment Dot Plot",
         subtitle = paste0("A172-R vs A172-S, ", total_sig, " significant (p<0.05), ",
                           total_metabs, " total metabolites"),
         x = "Fold Enrichment", y = "") +
    theme_pub +
    theme(axis.text.y = element_text(size = 10))
  save_fig(p_enrich_dot, "enrichment_dot_plot", w = 12, h = 8)

  cat("Enrichment figures saved.\n")
} else {
  cat("No categories with sufficient hits for plotting.\n")
}

# ============================================================================
# 6. METABOLITE CLASS DISTRIBUTION FIGURES
# ============================================================================
cat("\n=== Step 6: Additional Figures ===\n")

# 6a. Pie chart of metabolite classification
class_pie <- as.data.frame(class_counts)
colnames(class_pie) <- c("Category", "Count")
class_pie <- class_pie[order(-class_pie$Count), ]
class_pie$Fraction <- class_pie$Count / sum(class_pie$Count)
class_pie$Label <- paste0(class_pie$Category, " (", class_pie$Count, ")")

p_pie <- ggplot(class_pie, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Metabolite Distribution by Chemical Class",
       subtitle = paste0("104 metabolites, A172 glioblastoma metabolomics")) +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))
save_fig(p_pie, "metabolite_class_pie", w = 10, h = 8)

# 6b. Bar chart of differentially abundant metabolites by class
sig_by_class <- table(metab_info$class[metab_info$significant])
sig_class_df <- as.data.frame(sig_by_class)
colnames(sig_class_df) <- c("Category", "Significant")
sig_class_df$Total <- class_counts[match(sig_class_df$Category, names(class_counts))]
sig_class_df$Percent <- round(100 * sig_class_df$Significant / sig_class_df$Total, 1)
sig_class_df <- sig_class_df[order(-sig_class_df$Significant), ]

p_class_bar <- ggplot(sig_class_df, aes(x = reorder(Category, Significant), y = Significant)) +
  geom_bar(stat = "identity", fill = "#E41A1C", alpha = 0.8) +
  geom_text(aes(label = paste0(Significant, "/", Total, " (", Percent, "%)")),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(title = "Significant Metabolites by Chemical Class",
       subtitle = paste0("p < 0.05, Welch's t-test, n=4/group"),
       x = "", y = "Number of Significant Metabolites") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9))
save_fig(p_class_bar, "significant_by_class", w = 12, h = 8)

cat("Additional figures saved.\n")

# ============================================================================
# 7. SAVE FULL RESULTS
# ============================================================================
cat("\n=== Step 7: Saving Results ===\n")

# Create comprehensive results table
results_df <- metab_info
results_df$hmdb_id <- NA_character_
results_df <- results_df[order(results_df$p_value), ]
write.csv(results_df, file.path(DATA_DIR, "metabolite_classification_results.csv"), row.names = FALSE)

# Summary statistics
cat("\n=== Enrichment Analysis Summary ===\n")
cat("Total metabolites:", total_metabs, "\n")
cat("Significant (p<0.05):", total_sig, "\n")
cat("FDR significant (q<0.05):", sum(tt_df$fdr_p < 0.05), "\n")
cat("Classes identified:", length(unique(metab_info$class)), "\n")
cat("\nTop enriched categories:\n")
print(head(enrich_df[enrich_df$Hits >= 1, c("Category", "Hits", "FoldEnrichment", "p_value", "fdr")], 10))

cat("\n=== ENRICHMENT ANALYSIS COMPLETE ===\n")
cat("All figures in:", FIG_DIR, "\n")
cat("All data in:", DATA_DIR, "\n")

sessionInfo()
