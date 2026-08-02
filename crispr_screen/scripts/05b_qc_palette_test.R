#!/usr/bin/env Rscript
# Fig05 palette comparison — 4 candidate palettes, exported for visual selection.
# Saves palette1–4 variants to REVISED_DIR. The user selects one after inspection.
source("utils.R")
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

# ---- palette definitions -------------------------------------------------------
palettes <- list(
  palette1 = c(S0 = "#009E73", Spost = "#0072B2", R0 = "#D55E00", Rpost = "#CC79A7"),
  palette2 = c(S0 = "#4DBBD5", Spost = "#00A087", R0 = "#E64B35", Rpost = "#7E6148"),
  palette3 = c(S0 = "#8ECFC9", Spost = "#0072B2", R0 = "#F39B7F", Rpost = "#D55E00"),
  palette4 = c(S0 = "#1B9E77", Spost = "#386CB0", R0 = "#E7298A", Rpost = "#D95F02")
)

# ---- data loading (same as 05_qc_metrics.R) ------------------------------------
count_file <- file.path(DATA_DIR, "count_table.txt")
raw <- read_tsv(count_file, show_col_types = FALSE)

count_cols <- intersect(c("R0_1", "R0_2", "R0_3",
                          "Rpost_1", "Rpost_2", "Rpost_3", "Rpost_4",
                          "S0_1", "S0_2", "S0_3",
                          "Spost_1", "Spost_2", "Spost_3"), names(raw))

sample_cond <- sub("_[0-9]+$", "", count_cols)
names(sample_cond) <- count_cols
label_map <- sub("_", "-", count_cols)
names(label_map) <- count_cols

sample_order <- count_cols[order(match(sample_cond[count_cols], COND_LEVELS),
                                 as.integer(sub(".*_", "", count_cols)))]
label_order  <- unname(label_map[sample_order])

counts <- raw %>%
  select(all_of(count_cols)) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "count") %>%
  mutate(condition = factor(sample_cond[sample], levels = COND_LEVELS),
         sample_label = factor(label_map[sample], levels = label_order))

gini <- function(x) {
  x <- sort(x[x > 0]); n <- length(x)
  if (n < 2) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

# ---- build function ------------------------------------------------------------
build_qc_figure <- function(pal, pname) {
  pal_colors <- pal

  per_sample <- counts %>%
    group_by(sample_label, condition) %>%
    summarise(total    = sum(count, na.rm = TRUE) / 1e6,
              detected = sum(count > 0, na.rm = TRUE),
              .groups = "drop")

  pA <- ggplot(per_sample, aes(sample_label, total, fill = condition)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = pal_colors, guide = "none") +
    labs(x = NULL, y = "Toplam okuma (milyon)") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid.major.x = element_blank())

  pB <- ggplot(per_sample, aes(sample_label, detected, fill = condition)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = detected), vjust = -0.4, size = 2.6, family = BASE_FAMILY) +
    scale_y_log10(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = pal_colors, guide = "none") +
    labs(x = NULL, y = "Saptanan sgRNA") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid.major.x = element_blank())

  gini_df <- counts %>%
    group_by(sample_label, condition) %>%
    summarise(gini = gini(count), .groups = "drop")

  pC <- ggplot(gini_df, aes(sample_label, gini, fill = condition)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = pal_colors, guide = "none") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Gini indeksi") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid.major.x = element_blank())

  lorenz <- counts %>%
    filter(count > 0) %>%
    group_by(sample_label, condition) %>%
    arrange(count, .by_group = TRUE) %>%
    mutate(cum_reads = cumsum(count) / sum(count),
           cum_sgrna = row_number() / n()) %>%
    ungroup()

  pD <- ggplot(lorenz, aes(cum_sgrna * 100, cum_reads * 100,
                           colour = condition, group = sample_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    scale_colour_manual(values = pal_colors, labels = COND_LABELS, name = "Grup") +
    labs(x = "Kümülatif sgRNA (%)", y = "Kümülatif okuma (%)") +
    theme_thesis()

  combined <- (pA | pB) / (pC | pD) +
    plot_annotation(
      title    = "Kütüphane kalite kontrolü",
      subtitle = "30.197 sgRNA x 13 örnek; yüksek seyrekliğe dikkat edilmelidir.",
      tag_levels = "A"
    ) &
    tag_theme() &
    theme(legend.position = "bottom",
          plot.title    = element_text(margin = margin(l = 16, b = 2)),
          plot.subtitle = element_text(margin = margin(l = 16, b = 6)))

  save_figure(combined, paste0("Fig05_qc_metrics_", pname),
              width = 11, height = 8.5)
}

# ---- generate all 4 variants ---------------------------------------------------
for (nm in names(palettes)) {
  message("Generating ", nm, " ...")
  build_qc_figure(palettes[[nm]], nm)
}
message("All 4 palette variants saved to ", REVISED_DIR)
