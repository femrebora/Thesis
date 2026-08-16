# utils.R -- shared helpers for CRISPR screen figures
# PRJ25-131 | Cingoz Lab

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(grid)
})

# Deterministic label placement (ggrepel) and any stochastic geom across figures
set.seed(42)

# -- Typography ----------------------------------------------------------------
# Prefer an Arial-metric sans for a clean thesis look; fall back gracefully so
# the pipeline never errors on a machine without the exact font.
.pick_family <- function() {
  fams <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
  for (f in c("Liberation Sans", "Nimbus Sans", "Arial", "DejaVu Sans")) {
    if (f %in% fams) return(f)
  }
  "sans"
}
BASE_FAMILY <- .pick_family()

# -- Paths ---------------------------------------------------------------------
# All paths are anchored to the module directory via here::here(), so scripts run
# correctly regardless of the working directory. Inputs (count matrix, MAGeCK
# gene summaries) live in data/ and are git-ignored; see crispr_screen/README.md
# and crispr_screen/data/README.md.
.MODULE_DIR <- here::here("crispr_screen")
DATA_DIR    <- file.path(.MODULE_DIR, "data")
FIG_DIR     <- file.path(.MODULE_DIR, "figures")
REVISED_DIR <- FIG_DIR                        # canonical figure output (thesis-final versions)
TAB_DIR     <- file.path(.MODULE_DIR, "results")
GENE_DIR    <- file.path(DATA_DIR, "gene_lists")

# -- Significance criterion (configurable) -------------------------------------
# The thesis baseline calls a gene significant on its NOMINAL MAGeCK RRA p-value
# (p < 0.05, uncorrected). The whole figure set can be rebuilt under a
# Benjamini-Hochberg FDR acceptance criterion instead, without touching any
# figure script, by exporting:
#
#   THESIS_SIG_METRIC=fdr  THESIS_SIG_THRESH=0.05  THESIS_SIG_SUFFIX=_fdr05
#
# See scripts/run_crispr_analysis_fdr05.R. Outputs then carry SIG_SUFFIX so the
# baseline figures are never overwritten.
#
# NOTE: Several figure scripts still declare a local `pval_thresh <- 0.05` and
# pass that into loaders/plot code. At the thesis default this equals SIG_THRESH;
# it is not a fully centralised single constant across every script. See
# docs/methods_and_thresholds.md. Do not "fix" those locals without proving
# byte-identical PNG regeneration under the recorded environment.
SIG_METRIC <- tolower(Sys.getenv("THESIS_SIG_METRIC", "pvalue"))
if (!SIG_METRIC %in% c("pvalue", "fdr")) {
  stop("THESIS_SIG_METRIC must be 'pvalue' or 'fdr', got: ", SIG_METRIC)
}
SIG_THRESH <- as.numeric(Sys.getenv("THESIS_SIG_THRESH", "0.05"))
SIG_SUFFIX <- Sys.getenv("THESIS_SIG_SUFFIX", "")

# Turkish axis / caption token for whichever statistic is in force.
SIG_LABEL_TR <- if (SIG_METRIC == "fdr") "FDR" else "p-değeri"
SIG_SHORT_TR <- if (SIG_METRIC == "fdr") "FDR" else "p"

# -- Colour palette ------------------------------------------------------------
# Defined ONCE here and reused everywhere. Two independent palettes so a hue
# never carries two meanings across the figure set.
#
# (1) DIRECTION palette -- encodes the direction of an sgRNA representation
#     change only (Azalan / Artan / Anlamlı değil).
C_DEPLETED  <- "#2C7FB8"   # muted deep blue   -> Azalan (depleted)
C_ENRICHED  <- "#D95F02"   # muted orange      -> Artan  (enriched)
C_NS        <- "#BDBDBD"   # neutral grey      -> Anlamlı değil (not significant)
C_HIGHLIGHT <- "#F0E442"
C_NA        <- "#F2F2F2"   # very light grey   -> not measured / NA cells (heatmaps)

# (2) GROUP / CONDITION palette -- Sensitive = teal family, Resistant = purple
#     family (colourblind-safe, deliberately OFF the blue/orange direction hues).
#     Light shade = Day 0 baseline; dark shade = Tumor-derived (in vivo).
COND_COLORS <- c(
  S0    = "#009E73",   # Duyarlı  · Gün 0    (green — Okabe-Ito)
  Spost = "#0072B2",   # Duyarlı  · Tümör    (blue — Okabe-Ito)
  R0    = "#D55E00",   # Dirençli · Gün 0    (orange-red — Okabe-Ito)
  Rpost = "#CC79A7"    # Dirençli · Tümör    (reddish-purple — Okabe-Ito)
)

# Short academic-Turkish group labels for legends / axes.
COND_LABELS <- c(
  S0    = "Duyarlı Gün 0",
  Spost = "Duyarlı Tümör",
  R0    = "Dirençli Gün 0",
  Rpost = "Dirençli Tümör"
)
COND_LEVELS <- c("S0", "Spost", "R0", "Rpost")

# -- Comparison labels ---------------------------------------------------------
# Short academic-Turkish comparison labels (consistent across every figure).
# 'Tümör' = tumor-derived after in vivo growth (NO TRAIL treatment in this screen);
# 'Gün 0' = Day 0 baseline pool. Resistance is a cell-line property (A172 vs A172R).
COMPARISONS <- c(
  Rpost_vs_R0    = "Dirençli: Tümör / Gün 0",
  Spost_vs_Rpost = "Duyarlı / Dirençli (Tümör)",
  Spost_vs_S0    = "Duyarlı: Tümör / Gün 0"
)

# -- ggplot2 theme -------------------------------------------------------------
# Print-calibrated, font-embedded, no clipped text. See FIGURE_STANDARD.md.
theme_thesis <- function(base_size = 10, family = BASE_FAMILY) {
  theme_classic(base_size = base_size, base_family = family) +
    theme(
      text = element_text(family = family, color = "#1A1A1A"),
      panel.grid.major.y = element_line(color = "#EAEAEA", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "#333333", linewidth = 0.4),
      axis.ticks = element_line(color = "#333333", linewidth = 0.4),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "#333333", size = base_size - 1),
      axis.title = element_text(color = "#1A1A1A", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 2,
                                margin = margin(b = 4)),
      plot.subtitle = element_text(color = "#555555", size = base_size - 1,
                                   margin = margin(b = 6)),
      plot.title.position = "plot",
      plot.caption = element_text(color = "#777777", size = base_size - 2),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.9, "lines"),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 6, 8)
    )
}

# Patchwork tag styling (A, B, C ...) consistent across multi-panel figures.
tag_theme <- function(base_size = 10, family = BASE_FAMILY) {
  theme(plot.tag = element_text(face = "bold", size = base_size + 3,
                                family = family),
        plot.tag.position = c(0, 1))
}

# -- Data loader ---------------------------------------------------------------
load_gene_summary <- function(comparison, pval_thresh = 0.05) {
  fname <- paste0(comparison, ".gene_summary.tsv")
  path <- file.path(DATA_DIR, fname)
  if (!file.exists(path)) stop("Missing: ", path)

  df <- read_tsv(path, show_col_types = FALSE)

  # Raw MAGeCK statistics are always kept under *_raw / *_fdr names so any
  # script can reach either tail on either scale.
  df <- df %>%
    mutate(
      LFC       = `neg|lfc`,
      neg_p_raw = `neg|p-value`,          # depletion test (significant when LFC < 0)
      pos_p_raw = `pos|p-value`,          # enrichment test (significant when LFC > 0)
      neg_fdr   = `neg|fdr`,
      pos_fdr   = `pos|fdr`,
      pval_raw  = pmin(`neg|p-value`, `pos|p-value`),
      fdr       = pmin(`neg|fdr`, `pos|fdr`)
    )

  # neg_p / pos_p / pval carry the statistic CURRENTLY IN FORCE (see SIG_METRIC).
  # Downstream figure scripts read these three columns only, so switching the
  # acceptance criterion needs no change in any figure script.
  if (SIG_METRIC == "fdr") {
    df <- df %>% mutate(neg_p = neg_fdr, pos_p = pos_fdr, pval = fdr)
  } else {
    df <- df %>% mutate(neg_p = neg_p_raw, pos_p = pos_p_raw, pval = pval_raw)
  }

  df %>%
    mutate(
      log10p     = -log10(pmax(pval, 1e-100)),
      sig        = pval < pval_thresh,
      direction  = ifelse(LFC < 0, "Depleted", "Enriched"),
      comparison = COMPARISONS[comparison]
    )
}

# -- Save helper -- PDF (cairo, fonts embedded) + 600-dpi PNG ------------------
# Revised thesis outputs go to REVISED_DIR with a "_revised" suffix so the
# ORIGINAL figures in FIG_DIR are never overwritten.
save_figure <- function(plot, name, width = 7, height = 5, dpi = 600,
                        outdir = REVISED_DIR, suffix = "") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stem <- paste0(name, suffix, SIG_SUFFIX)

  ggsave(file.path(outdir, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = dpi, bg = "white",
         type = if (capabilities("cairo")) "cairo" else NULL)

  message("Saved: ", stem, " (.pdf + .png)  [", width, "x", height, " in]")
}

# Save a base-graphics / grid figure (e.g. ComplexHeatmap) to PDF + PNG with the
# same conventions as save_figure(). `draw_fn` is a function that renders the plot.
save_grid_figure <- function(draw_fn, name, width = 7, height = 5, dpi = 600,
                             outdir = REVISED_DIR, suffix = "") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stem <- paste0(name, suffix, SIG_SUFFIX)

  cairo_pdf(file.path(outdir, paste0(stem, ".pdf")), width = width, height = height)
  draw_fn(); dev.off()

  png(file.path(outdir, paste0(stem, ".png")), width = width, height = height,
      units = "in", res = dpi, bg = "white",
      type = if (capabilities("cairo")) "cairo" else "Xlib")
  draw_fn(); dev.off()

  message("Saved: ", stem, " (.pdf + .png)  [", width, "x", height, " in]")
}

# -- Figure helpers ------------------------------------------------------------
# Wrap long category / facet labels so nothing is clipped.
wrap_labels <- function(x, width = 24) str_wrap(x, width = width)

# Ordered Turkish direction factor levels (single source of truth).
DIR_LEVELS <- c("Azalan", "Artan", "Anlamlı değil")

# Sign-based classification (significance = overall min-p). Kept for figures that
# only need a coarse direction split.
classify_direction <- function(lfc, sig) {
  factor(ifelse(!sig, "Anlamlı değil",
                ifelse(lfc < 0, "Azalan", "Artan")),
         levels = DIR_LEVELS)
}

# Direction-CONCORDANT classification (statistically correct MAGeCK reading):
#  - Azalan : depletion test significant (neg_p < thresh) AND LFC < 0
#  - Artan  : enrichment test significant (pos_p < thresh) AND LFC > 0
#  - else   : Anlamlı değil  (includes discordant genes significant on the
#             opposite tail to their LFC sign -- e.g. TECR/TKTL1/APOC1).
# Recolours discordant points WITHOUT removing any rows from the data.
classify_direction_concordant <- function(lfc, neg_p, pos_p, thresh = 0.05) {
  dep <- !is.na(neg_p) & neg_p < thresh & lfc < 0
  enr <- !is.na(pos_p) & pos_p < thresh & lfc > 0
  factor(ifelse(dep, "Azalan", ifelse(enr, "Artan", "Anlamlı değil")),
         levels = DIR_LEVELS)
}

# Named colour scale matching the Turkish direction levels.
DIRECTION_COLORS <- c(
  "Azalan"        = C_DEPLETED,
  "Artan"         = C_ENRICHED,
  "Anlamlı değil" = C_NS
)

# -- Metabolic library universe (background for enrichment) --------------------
# The statistically correct enrichment background for a focused screen is the
# set of genes the library actually targets, NOT the whole genome.
load_library_universe <- function() {
  ct <- read_tsv(file.path(DATA_DIR, "count_table.txt"),
                 show_col_types = FALSE)
  genes <- unique(ct$Gene)
  genes <- genes[!grepl("^(INTERGENIC|CONTROL|NonTargeting|NTC|SafeHarbor|Safe)",
                        genes, ignore.case = TRUE)]
  genes[!is.na(genes) & nzchar(genes)]
}

# -- Gene annotation helpers ---------------------------------------------------
KNOWN_DEPLETED <- c("CYP2E1", "ASL", "CYB5A", "MTF1", "NAT2", "GBA3",
                    "TRPM7", "ATP6V0A4", "TPCN1", "HSD11B1", "PFKFB1")
KNOWN_ENRICHED <- c("PLA2G4E", "SLC1A1", "ARSD", "GRIN1", "RRM2", "SLC6A3",
                    "SLC44A1", "DHRSX", "ALOXE3", "PTGS1", "SLC25A15")

FUNC_ANNOT <- c(
  ASL      = "Argininosuccinate lyase - urea cycle, arginine biosynthesis",
  PLA2G4E  = "Cytosolic PLA2 - arachidonate release, inflammation",
  SLC1A1   = "EAAT3 glutamate transporter, GSH synthesis",
  GBA3     = "Cytosolic beta-glucosidase - glycolipid catabolism",
  RRM2     = "Ribonucleotide reductase M2 - dNTP synthesis, DNA repair",
  ATP6V0A4 = "V-ATPase subunit - lysosomal acidification",
  ARSD     = "Arylsulfatase D - lysosomal sulfatase",
  TPCN1    = "Two-pore channel 1 - lysosomal Ca2+ release",
  GRIN1    = "NMDA receptor GluN1 - Ca2+ signalling",
  HSD11B1  = "11beta-HSD1 - glucocorticoid activation",
  PFKFB1   = "6-phosphofructo-2-kinase - glycolysis regulation",
  TRPM7    = "TRPM7 channel - Mg2+/Ca2+ homeostasis",
  ALOXE3   = "Epidermal lipoxygenase-3 - ceramide signalling",
  SLC25A15 = "Mitochondrial ornithine carrier - urea cycle",
  SLC6A3   = "Dopamine transporter (DAT)",
  CYB5A    = "Cytochrome b5A - electron transfer",
  DHRSX    = "Dehydrogenase/reductase X",
  SLC44A1  = "Choline transporter-like protein 1",
  NAT2     = "N-acetyltransferase 2 - drug/xenobiotic acetylation",
  CYP2E1   = "CYP450 2E1 - xenobiotic metabolism, ROS generation",
  ALOX15   = "Arachidonate 15-lipoxygenase - lipid peroxidation, ferroptosis",
  SLC25A20 = "Carnitine/acylcarnitine carrier - mitochondrial fatty-acid oxidation"
)

get_annotation <- function(gene) {
  if (gene %in% names(FUNC_ANNOT)) FUNC_ANNOT[gene] else "-"
}

# Turkish function annotations (faithful translations of FUNC_ANNOT) for the
# Turkish-language thesis tables. Biological descriptions only -- no new claims.
FUNC_ANNOT_TR <- c(
  ASL      = "Üre döngüsü; arginin biyosentezi",
  PLA2G4E  = "Fosfolipaz A2 aktivitesi; araşidonat salınımı ve lipid aracıları",
  SLC1A1   = "Glutamat taşıyıcılığı; sistein alımı ve glutatyon dengesi",
  GBA3     = "Sitozolik beta-glukozidaz; glikozid ve glikolipid metabolizması",
  RRM2     = "dNTP sentezi; DNA replikasyonu ve onarımı",
  ATP6V0A4 = "V-ATPaz alt birimi; organel asidifikasyonu",
  ARSD     = "Sülfataz aktivitesi; lizozomal sülfat metabolizması",
  TPCN1    = "Endolizozomal iyon kanalı; Ca2+ sinyallemesi",
  GRIN1    = "NMDA reseptör alt birimi; glutamat aracılı Ca2+ sinyali",
  HSD11B1  = "11beta-HSD1 – glukokortikoid aktivasyonu",
  PFKFB1   = "6-fosfofrukto-2-kinaz – glikoliz düzenlenmesi",
  TRPM7    = "TRPM7 kanalı – Mg2+/Ca2+ homeostazı",
  ALOXE3   = "Epidermal lipoksigenaz-3 – seramid sinyalizasyonu",
  SLC25A15 = "Mitokondriyal ornitin taşıyıcısı – üre döngüsü",
  SLC6A3   = "Dopamin taşıyıcısı (DAT)",
  CYB5A    = "Sitokrom b5A – elektron transferi",
  DHRSX    = "Dehidrogenaz/redüktaz X",
  SLC44A1  = "Kolin taşıyıcı benzeri protein 1",
  NAT2     = "Ksenobiyotik ve ilaç asetilasyonu",
  CYP2E1   = "CYP450 2E1 – ksenobiyotik metabolizması, ROS üretimi",
  ALOX15   = "Lipoksigenaz aktivitesi; lipid peroksidasyonu ve eikosanoid metabolizması",
  SLC25A20 = "Karnitin-asilkarnitin taşınımı; mitokondriyal yağ asidi oksidasyonu"
)

get_annotation_tr <- function(gene) {
  if (gene %in% names(FUNC_ANNOT_TR)) FUNC_ANNOT_TR[gene] else "Açıklama yok"
}
