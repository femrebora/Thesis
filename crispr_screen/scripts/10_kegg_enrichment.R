#!/usr/bin/env Rscript
# Fig10 -- KEGG pathway enrichment (v3): metabolic-focused horizontal barplot.
# TÜSEB-inspired style: warm red bars, no gene counts, no main title.
# Non-metabolic terms (addiction, neurodegenerative, neuronal, infectious) excluded.
# Uses pre-computed enrichment tables (no API calls).
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

# ---- load data ----------------------------------------------------------------
kegg_all   <- read_csv(file.path(TAB_DIR, "kegg_all_25.csv"), show_col_types = FALSE)
kegg_rpost <- read_csv(file.path(TAB_DIR, "kegg_Rpost_18.csv"), show_col_types = FALSE)
kegg_all$GeneSet   <- "all_25"
kegg_rpost$GeneSet <- "Rpost_18"
kegg_df <- bind_rows(kegg_all, kegg_rpost)

# ---- filter: exclude non-metabolic terms ---------------------------------------
NON_METABOLIC_PATTERNS <- c(
  "addiction", "cocaine", "amphetamine", "alcoholism", "nicotine", "morphine",
  "Parkinson", "Alzheimer", "Huntington", "prion", "spinocerebellar",
  "neurodegeneration", "synaptic vesicle", "glutamatergic synapse",
  "serotonergic synapse", "dopaminergic synapse", "GABAergic synapse",
  "cholinergic synapse", "long-term potentiation", "long-term depression",
  "circadian", "neuroactive ligand", "neurotrophin", "axon guidance",
  "viral", "hepatitis", "influenza", "tuberculosis", "salmonella",
  "shigella", "legionella", "yersinia", "pertussis", "malaria",
  "leishmania", "toxoplasma", "amoebiasis", "cholera", "pathogenic",
  "HTLV", "HPV", "EBV", "CMV", "Kaposi", "measles", "influenza",
  "glioma", "leukemia", "breast cancer", "prostate cancer",
  "pancreatic cancer", "colorectal", "gastric cancer", "hepatocellular",
  "renal cell", "thyroid cancer", "bladder cancer", "endometrial",
  "basal cell", "small cell lung", "non-small cell lung",
  "melanoma", "chronic myeloid", "acute myeloid",
  "necroptosis", "ferroptosis", "cellular senescence",
  "Rap1 signaling", "Ras signaling", "MAPK signaling", "cAMP signaling",
  "calcium signaling", "AMPK signaling", "VEGF signaling",
  "phospholipase D signaling", "NOD-like receptor", "Fc epsilon",
  "Fc gamma", "platelet activation", "phagocytosis", "efferocytosis",
  "GnRH", "oxytocin", "glucagon", "insulin secretion", "insulin signaling",
  "mineral absorption", "protein digestion", "ovarian steroidogenesis",
  "collecting duct acid", "thermogenesis", "inflammatory mediator",
  "vascular smooth muscle", "cornified envelope",
  "choline metabolism in cancer", "central carbon metabolism in cancer",
  "cGMP-PKG", "apelin", "cortisol", "cushing", "estrogen",
  "parathyroid", "renin", "aldosterone", "thyroid hormone",
  "relaxin", "prolactin", "growth hormone", "adipocytokine",
  "type I diabetes", "type II diabetes", "maturity onset diabetes",
  "AGE-RAGE", "HIF-1", "FoxO", "PI3K-Akt", "JAK-STAT",
  "NF-kappa B", "TNF", "Toll-like", "RIG-I-like", "C-type lectin",
  "Th1 and Th2", "Th17", "T cell receptor", "B cell receptor",
  "chemokine", "complement and coagulation", "neutrophil extracellular",
  "antigen processing", "natural killer", "IL-17",
  "osteoclast", "rheumatoid arthritis", "graft-versus-host",
  "allograft rejection", "autoimmune thyroid", "systemic lupus",
  "asthma", "inflammatory bowel", "malaria", "African trypanosomiasis",
  "Chagas", "leishmaniasis", "toxoplasmosis", "amoebiasis",
  "Staphylococcus", "vibrio cholerae", "Helicobacter pylori",
  "Epstein-Barr", "human cytomegalovirus", "human papillomavirus",
  "Kaposi sarcoma", "hepatitis C", "hepatitis B", "measles",
  "influenza A", "COVID-19", "coronavirus",
  "morphine addiction", "GABAergic", "retrograde endocannabinoid",
  "dilated cardiomyopathy", "hypertrophic cardiomyopathy",
  "arrhythmogenic right ventricular", "cardiac muscle contraction",
  "adrenergic signaling", "vascular smooth muscle contraction",
  "gap junction", "tight junction", "adherens junction", "focal adhesion",
  "ECM-receptor", "cell adhesion molecules", "leukocyte transendothelial",
  "regulation of actin cytoskeleton", "apoptosis", "p53",
  "cell cycle", "oocyte meiosis", "progesterone-mediated",
  "cellular senescence", "lysosome", "peroxisome", "phagosome",
  "protein processing in endoplasmic reticulum", "SNARE interactions",
  "ubiquitin mediated", "spliceosome", "mRNA surveillance",
  "RNA degradation", "RNA polymerase", "DNA replication",
  "base excision", "nucleotide excision", "mismatch repair",
  "homologous recombination", "fanconi anemia", "non-homologous end-joining",
  "basal transcription factors", "ribosome", "ribosome biogenesis",
  "protein export", "proteasome", "RNA transport",
  "sphingolipid signaling", "sphingolipid metabolism"
)
# NOTE: "sphingolipid metabolism" is excluded above because the KEGG term
# "Sphingolipid metabolism" appears in non-metabolic context. The metabolic
# sphingolipid terms appear in Reactome instead.

kegg_metabolic <- kegg_df %>%
  filter(!grepl(paste(NON_METABOLIC_PATTERNS, collapse = "|"),
                Description, ignore.case = TRUE))

# ---- panel labels --------------------------------------------------------------
GENESET_TR <- c(
  all_25   = "Tüm genler",
  Rpost_18 = "Dirençli: Tümör / Gün 0"
)
kegg_metabolic$GeneSetTR <- factor(GENESET_TR[kegg_metabolic$GeneSet],
                                    levels = unname(GENESET_TR))

# ---- top 10 per panel -----------------------------------------------------------
kegg_top <- kegg_metabolic %>%
  group_by(GeneSetTR) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    log10padj   = -log10(p.adjust),
    Description = str_wrap(Description, width = 45)
  )

kegg_top <- kegg_top %>%
  group_by(GeneSetTR) %>%
  mutate(Description = forcats::fct_reorder(Description, log10padj)) %>%
  ungroup()

# ---- plot ----------------------------------------------------------------------
p <- ggplot(kegg_top, aes(x = log10padj, y = Description, fill = log10padj)) +
  geom_col(width = 0.7, alpha = 0.9) +
  scale_fill_gradientn(
    colours = c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
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

save_figure(p, "Fig10_kegg_enrichment_v3", suffix = "",
            width = 8, height = max(6, length(unique(kegg_top$Description)) * 0.33 + 2))
