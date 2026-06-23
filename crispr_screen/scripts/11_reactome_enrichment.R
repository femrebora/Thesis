#!/usr/bin/env Rscript
# Fig11 -- Reactome pathway enrichment (v3): metabolic-focused horizontal barplot.
# TÜSEB-inspired style: warm red bars, no gene counts, no main title.
# Non-metabolic terms (neuronal, signaling, disease) excluded.
# Uses pre-computed enrichment tables (no API calls).
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

# ---- load data ----------------------------------------------------------------
reac_all   <- read_csv(file.path(TAB_DIR, "reactome_all_25.csv"), show_col_types = FALSE)
reac_rpost <- read_csv(file.path(TAB_DIR, "reactome_Rpost_18.csv"), show_col_types = FALSE)
reac_all$GeneSet   <- "all_25"
reac_rpost$GeneSet <- "Rpost_18"
reac_df <- bind_rows(reac_all, reac_rpost)

# ---- filter: exclude non-metabolic terms ---------------------------------------
NON_METABOLIC_PATTERNS <- c(
  "NMDA", "glutamate", "Neurotransmitter", "Synaptic", "synapse",
  "Chemical Synapses", "Neuronal System", "Long-term potentiation",
  "Neurexins", "neuroligins", "Postsynaptic", "Post NMDA",
  "Unblocking of NMDA", "Activation of NMDA", "CREB1",
  "Ras activation upon Ca2+", "Assembly and cell surface presentation of NMDA",
  "Negative regulation of NMDA", "Protein-protein interactions at synapses",
  "Ion channel transport", "Stimuli-sensing channels", "TRP channels",
  "SLC-mediated transmembrane transport", "SLC transporter disorders",
  "Na+/Cl- dependent neurotransmitter",
  "Transmission across Chemical Synapses",
  "Insulin receptor recycling", "Signaling by Insulin receptor",
  "EPH-Ephrin", "EPHB-mediated",
  "MAPK", "MAPK1", "MAPK3", "RAF", "Signaling by RAF",
  "G1/S", "G1/S-Specific Transcription", "Transcriptional Regulation by E2F6",
  "ROS and RNS production in phagocytes",
  "Transferrin endocytosis", "Iron uptake and transport",
  "Plasma lipoprotein", "lipoprotein assembly", "lipoprotein clearance",
  "NR1H2", "NR1H3", "Signaling by Nuclear Receptors",
  "Signaling by Interleukins", "IL-4", "IL-13",
  "Protein localization", "Insertion of tail-anchored",
  "Diseases", "disease", "disorders", "Disorders",
  "Defects", "defects", "Mutant", "mutant",
  "Viral", "viral", "Infection", "infection",
  "SARS-CoV", "COVID", "HIV", "Influenza",
  "Oncogenic", "Tumor suppressor", "Cancer"
)

reac_metabolic <- reac_df %>%
  filter(!grepl(paste(NON_METABOLIC_PATTERNS, collapse = "|"),
                Description, ignore.case = TRUE))

# ---- panel labels --------------------------------------------------------------
GENESET_TR <- c(
  all_25   = "Tüm genler",
  Rpost_18 = "Dirençli: Tümör / Gün 0"
)
reac_metabolic$GeneSetTR <- factor(GENESET_TR[reac_metabolic$GeneSet],
                                    levels = unname(GENESET_TR))

# ---- top 10 per panel -----------------------------------------------------------
reac_top <- reac_metabolic %>%
  group_by(GeneSetTR) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    log10padj   = -log10(p.adjust),
    Description = str_wrap(Description, width = 48)
  )

reac_top <- reac_top %>%
  group_by(GeneSetTR) %>%
  mutate(Description = forcats::fct_reorder(Description, log10padj)) %>%
  ungroup()

# ---- plot ----------------------------------------------------------------------
p <- ggplot(reac_top, aes(x = log10padj, y = Description, fill = log10padj)) +
  geom_col(width = 0.7, alpha = 0.9) +
  scale_fill_gradientn(
    colours = c("#FEE0D2", "#FC9272", "#DE2D26", "#A50F15"),
    guide = "none"
  ) +
  facet_wrap(~ GeneSetTR, ncol = 1, scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = expression(-Log[10]~"(düzeltilmiş p-değeri)"),
    y = NULL
  ) +
  theme_thesis() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_blank()
  )

save_figure(p, "Fig11_reactome_enrichment_v3", suffix = "",
            width = 8, height = max(6, length(unique(reac_top$Description)) * 0.33 + 2))
