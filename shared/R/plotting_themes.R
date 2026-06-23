# shared/R/plotting_themes.R — Cross-module publication-quality ggplot2 themes
# Reusable across tcga_gbm, metabolomics, and crispr_screen modules

suppressPackageStartupMessages({
  library(ggplot2)
})

# ============================================================================
# Wong (2011) colorblind-safe 7-color palette
# ============================================================================
wong_palette <- c(
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#009E73",  # green
  "#F0E442",  # yellow
  "#CC79A7",  # reddish purple
  "#56B4E9",  # sky blue
  "#E69F00"   # orange
)

# ============================================================================
# Okabe-Ito colorblind-safe 8-color palette
# ============================================================================
okabe_ito_palette <- c(
  "#000000",  # black
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7"   # reddish purple
)

# ============================================================================
# Publication journal theme (colorblind-safe, clean)
# ============================================================================
theme_thesis <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      text             = element_text(family = "sans", color = "#333333"),
      plot.title       = element_text(face = "bold", size = base_size + 3,
                                      hjust = 0, margin = margin(b = 8)),
      plot.subtitle    = element_text(size = base_size, hjust = 0,
                                      color = "#555555", margin = margin(b = 10)),
      axis.title       = element_text(size = base_size + 1),
      axis.title.x     = element_text(margin = margin(t = 6)),
      axis.title.y     = element_text(margin = margin(r = 6)),
      axis.text        = element_text(size = base_size - 1, color = "#444444"),
      axis.ticks       = element_line(color = "#cccccc"),
      legend.position  = "bottom",
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.box.spacing = unit(0, "mm"),
      legend.margin    = margin(t = 0),
      panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "#cccccc", fill = NA),
      strip.background = element_rect(fill = "#f5f5f5", color = "#cccccc"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = margin(12, 12, 12, 12)
    )
}

# Lightweight theme (for smaller panels, insets)
theme_thesis_light <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.3, color = "grey90"),
      axis.title        = element_text(size = base_size + 2, face = "bold"),
      axis.text         = element_text(size = base_size, color = "black"),
      legend.title      = element_text(face = "bold"),
      legend.position   = "bottom"
    )
}

message("[shared/plotting_themes] Loaded wong_palette, okabe_ito_palette, ",
        "theme_thesis, theme_thesis_light")
