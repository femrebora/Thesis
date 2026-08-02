# analysis/config.R
# Single source of truth for all analysis modules

# --- Gene Lists ---
GENES_DEPLETED <- c("ASL", "GBA3", "ATP6V0A4", "TPCN1")
GENES_ENRICHED <- c("PLA2G4E", "SLC1A1", "RRM2", "ARSD", "GRIN1")
GENES_OTHER    <- c("NAT2", "SLC25A20", "ALOX15")
ALL_GENES      <- c(GENES_DEPLETED, GENES_ENRICHED, GENES_OTHER)

# --- Classifications ---
GENE_CLASS <- c(
  setNames(rep("Depleted", length(GENES_DEPLETED)), GENES_DEPLETED),
  setNames(rep("Enriched", length(GENES_ENRICHED)), GENES_ENRICHED),
  setNames(rep("Other", length(GENES_OTHER)), GENES_OTHER)
)

# --- LFC values from thesis Table 4.2 ---
GENE_LFC_RESISTANT <- c(
  "PLA2G4E" = 12.3, "SLC1A1" = 11.4, "RRM2" = 10.9, "ARSD" = 10.4,
  "GRIN1" = 10.3, "NAT2" = 3.2, "SLC25A20" = -0.0, "TPCN1" = -6.7,
  "ATP6V0A4" = -7.3, "GBA3" = -8.8, "ASL" = -11.5, "ALOX15" = NA
)

# --- Paths ---
RESULTS_DIR <- "../results"
CACHE_DIR <- "GDCdata"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Module 01 Output ---
OUTPUT_DIR_01 <- file.path(RESULTS_DIR, "01_expression_survival")
FIGURES_DIR_01 <- file.path(OUTPUT_DIR_01, "figures")
TABLES_DIR_01  <- file.path(OUTPUT_DIR_01, "tables")

# --- TCGA Project ---
TCGA_PROJECT <- "TCGA-GBM"

# --- Caching ---
USE_CACHE <- TRUE

# --- Analysis thresholds ---
MIN_PATIENTS_KM       <- 20   # Minimum matched patients for survival analysis
MIN_NORMAL_SAMPLES    <- 5    # Minimum normal samples for tumor-normal comparison
MIN_GENES_HEATMAP     <- 8    # Minimum genes found for expression heatmap
MIN_EVENTS_PER_GROUP  <- 3    # Minimum events per KM group

# --- Output formats ---
FIGURE_FORMATS <- c("pdf", "png", "tiff")
FIGURE_DPI_PNG  <- 600
FIGURE_DPI_TIFF <- 600
TABLE_FORMATS  <- c("csv", "tsv", "xlsx")

message("[config] Loaded ", length(ALL_GENES), " genes (Depleted=", length(GENES_DEPLETED),
        " Enriched=", length(GENES_ENRICHED), " Other=", length(GENES_OTHER), ")")
