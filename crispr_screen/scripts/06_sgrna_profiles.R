#!/usr/bin/env Rscript
# Fig06 -- sgRNA count profiles of selected hit genes (rebuilt from scratch, v3).
# Design: one gene per facet, jittered individual sgRNA points + group mean line.
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(patchwork)
})

# ---- data ----------------------------------------------------------------------
HIT_GENES <- c("ASL", "PLA2G4E", "SLC1A1", "RRM2", "CYP2E1", "GRIN1")

count_file <- file.path(DATA_DIR, "count_table.txt")
raw <- read_tsv(count_file, show_col_types = FALSE)

count_cols <- intersect(c("R0_1", "R0_2", "R0_3",
                          "Rpost_1", "Rpost_2", "Rpost_3", "Rpost_4",
                          "S0_1", "S0_2", "S0_3",
                          "Spost_1", "Spost_2", "Spost_3"), names(raw))

sample_cond <- sub("_[0-9]+$", "", count_cols)
names(sample_cond) <- count_cols

# Short x-axis labels with controlled line breaks
COND_SHORT <- c(
  S0    = "Duyarlı\nGün 0",
  Spost = "Duyarlı\nTümör",
  R0    = "Dirençli\nGün 0",
  Rpost = "Dirençli\nTümör"
)

# ---- build long-format data for all 6 genes ------------------------------------
profiles <- raw %>%
  select(sgRNA, Gene, all_of(count_cols)) %>%
  filter(Gene %in% HIT_GENES) %>%
  pivot_longer(all_of(count_cols), names_to = "sample", values_to = "count") %>%
  mutate(
    condition = factor(sample_cond[sample], levels = COND_LEVELS),
    cond_label = factor(COND_SHORT[as.character(sample_cond[sample])],
                        levels = unname(COND_SHORT)),
    log_count = log10(count + 1),
    Gene = factor(Gene, levels = HIT_GENES)
  )

# Per-gene sgRNA count for facet subtitle
gene_labels <- profiles %>%
  group_by(Gene) %>%
  summarise(n_sg = n_distinct(sgRNA), .groups = "drop") %>%
  mutate(label = sprintf("%s  (%d sgRNA)", Gene, n_sg))
names(gene_labels$label) <- gene_labels$Gene

# ---- plot ----------------------------------------------------------------------
p <- ggplot(profiles, aes(x = cond_label, y = log_count, colour = condition)) +
  # Individual sgRNAs: subtle jittered points
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 42),
             size = 1.5, alpha = 0.45) +
  # Group mean: bold connecting line + white-filled coloured-border point
  stat_summary(aes(group = condition), fun = mean, geom = "line",
               linewidth = 1.0) +
  stat_summary(aes(group = condition), fun = mean, geom = "point",
               shape = 21, size = 3.0, fill = "white", stroke = 0.8) +
  scale_colour_manual(values = COND_COLORS, labels = COND_LABELS,
                      name = "Grup") +
  facet_wrap(~ Gene, ncol = 3, scales = "free_y") +
  labs(
    title    = "Seçilmiş genlerin sgRNA profilleri",
    y        = expression(Log[10]~"(okuma + 1)"),
    x        = NULL
  ) +
  theme_thesis() +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.major = element_blank()
  )

save_figure(p, "Fig06_sgrna_profiles_rebuilt", suffix = "", width = 11, height = 7.5)
