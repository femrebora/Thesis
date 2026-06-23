#!/usr/bin/env Rscript
# Fig12 -- GO enrichment (v3): three separate figures per ontology (BP, CC, MF).
# TÜSEB-inspired style: warm red bars, no gene counts, no main title.
# Metabolic/transporter-focused filtering per ontology.
# Uses pre-computed enrichment tables (no API calls).
source(here::here("crispr_screen", "R", "utils.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

# ---- load data ----------------------------------------------------------------
go_all   <- read_csv(file.path(TAB_DIR, "go_all_25.csv"), show_col_types = FALSE)
go_rpost <- read_csv(file.path(TAB_DIR, "go_Rpost_18.csv"), show_col_types = FALSE)
go_all$GeneSet   <- "all_25"
go_rpost$GeneSet <- "Rpost_18"
go_df <- bind_rows(go_all, go_rpost)

# ---- panel labels --------------------------------------------------------------
GENESET_TR <- c(
  all_25   = "Tüm genler",
  Rpost_18 = "Dirençli: Tümör / Gün 0"
)
go_df$GeneSetTR <- factor(GENESET_TR[go_df$GeneSet],
                           levels = unname(GENESET_TR))

# ---- filter lists per ontology -------------------------------------------------
# BP: exclude neuronal, addiction, behavioral, non-metabolic signaling
BP_EXCLUDE <- c(
  "dopamine", "catecholamine", "cocaine", "amphetamine", "morphine",
  "ethanol", "alcohol", "nicotine",
  "synaptic", "synapse", "neurotransmitter", "glutamate",
  "learning", "cognition", "behavior", "behaviour", "locomotory",
  "startle", "grooming", "prepulse", "reflex",
  "viral", "apoptotic process", "motor neuron",
  "protein tetramerization", "protein heterotetramerization",
  "protein heterooligomerization", "protein oligomerization",
  "endocytosis", "recycling", "regulation of endocytic",
  "response to mercury", "response to alkaloid",
  "response to organophosphorus", "response to cAMP",
  "response to ethanol", "response to drug",
  "sensory perception", "detection of", "visual",
  "calcium ion import", "calcium ion transport",
  "regulation of postsynaptic", "positive regulation of synaptic",
  "negative regulation of synaptic",
  "neuron projection", "axon", "dendrite", "dendritic",
  "cellular response to hormone", "response to hormone",
  "cellular response to endogenous", "response to endogenous",
  "cellular response to organic cyclic", "response to organic cyclic",
  "cellular response to nitrogen", "response to nitrogen compound",
  "cellular response to oxygen", "response to oxygen",
  "regulation of membrane potential", "regulation of ion",
  "regulation of transport", "positive regulation of transport",
  "negative regulation of transport", "regulation of secretion",
  "regulation of signaling", "positive regulation of signaling",
  "negative regulation of signaling",
  "signal transduction", "cell surface receptor signaling",
  "intracellular signal transduction", "second-messenger-mediated",
  "G protein-coupled", "adenylate cyclase", "phospholipase C",
  "protein kinase", "phosphatidylinositol", "MAPK cascade",
  "regulation of cell", "positive regulation of cell",
  "negative regulation of cell", "regulation of growth",
  "regulation of development", "regulation of differentiation",
  "regulation of proliferation", "regulation of apoptosis",
  "immune", "immune system", "inflammatory", "defense",
  "wound healing", "blood", "coagulation", "hemostasis",
  "circulatory", "heart", "cardiac", "vasculature",
  "kidney", "renal", "lung", "respiratory",
  "bone", "skeletal", "cartilage", "connective tissue",
  "epidermis", "skin", "hair",
  "digestive", "gut", "intestinal",
  "reproductive", "gonad", "ovarian", "testicular",
  "embryo", "embryonic", "morphogenesis", "pattern specification",
  "anatomical structure", "tissue", "organ"
)

# MF: exclude channel, receptor, non-metabolic binding
MF_EXCLUDE <- c(
  "channel activity", "voltage-gated", "ligand-gated",
  "calcium channel", "sodium channel", "potassium channel",
  "chloride channel", "cation channel", "anion channel",
  "neurotransmitter transmembrane transporter",
  "neurotransmitter receptor", "glutamate receptor",
  "NMDA", "AMPA", "GABA", "dopamine binding",
  "dopamine", "catecholamine", "serotonin",
  "G protein-coupled", "cytokine receptor",
  "growth factor receptor", "receptor ligand",
  "receptor signaling", "hormone receptor",
  "Hsp90", "Hsp70", "heat shock protein",
  "syntaxin", "SNARE", "syntaxin binding",
  "phosphatidylinositol bisphosphate binding",
  "phosphatidylinositol-3,5-bisphosphate",
  "phosphatidylinositol phosphate binding",
  "quaternary ammonium group binding",
  "structural constituent", "structural molecule",
  "extracellular matrix", "cytoskeletal",
  "motor activity", "actin binding", "actin filament",
  "microtubule", "tubulin binding",
  "DNA binding", "RNA binding", "transcription",
  "translation", "nucleic acid",
  "nucleotide binding", "ribonucleotide binding",
  "protein binding", "kinase activity",
  "protein kinase", "protein phosphatase",
  "enzyme activator", "enzyme inhibitor",
  "signaling receptor", "receptor activity",
  "molecular transducer", "molecular adaptor",
  "protein-containing complex", "protein dimerization",
  "protein homodimerization", "protein heterodimerization",
  "identical protein binding"
)

# CC: exclude neuronal compartments, keep metabolic/transport-related
CC_EXCLUDE <- c(
  "axon", "synaptic", "synapse", "neuron projection",
  "dendritic", "dendrite", "postsynaptic", "presynaptic",
  "glutamatergic", "GABAergic", "dopaminergic", "cholinergic",
  "NMDA", "AMPA", "voltage-gated", "sodium channel complex",
  "calcium channel complex", "potassium channel complex",
  "myelin", "node of Ranvier",
  "nuclear", "nucleus", "nucleoplasm", "nucleolus", "chromosome",
  "chromatin", "nuclear envelope", "nuclear pore",
  "ribosome", "ribosomal", "cytosolic ribosome",
  "spliceosomal", "spliceosome", "catalytic step 2",
  "U2-type", "U4/U6", "U5", "SMN",
  "transcription", "transcriptional", "mediator",
  "RNA polymerase", "PcG", "PRC1", "histone",
  "proteasome", "proteasomal",
  "centrosome", "centriolar", "centriole", "spindle",
  "microtubule", "cytoskeleton", "actin", "myosin",
  "intermediate filament", "keratin",
  "extracellular", "extracellular matrix", "collagen",
  "basement membrane", "laminin", "fibronectin",
  "focal adhesion", "cell-cell junction", "tight junction",
  "adherens junction", "desmosome", "gap junction",
  "apical", "basolateral", "lateral plasma",
  "cell surface", "external side",
  "brush border", "microvillus", "stereocilium",
  "cilium", "ciliary", "axoneme", "flagellum",
  "sperm", "acrosomal", "male",
  "blood", "hemoglobin", "platelet", "serum",
  "immunological", "immunoglobulin", "MHC",
  "T cell", "B cell", "TCR", "BCR",
  "secretory", "secretory granule", "zymogen",
  "vesicle", "synaptic vesicle", "transport vesicle",
  "COPI", "COPII", "clathrin", "caveola",
  "melanosome", "pigment",
  "lipid droplet", "droplet",
  "contractile", "sarcomere", "myofibril", "Z disc", "I band",
  "sarcoplasmic", "T-tubule"
)

# ---- build one figure per ontology ---------------------------------------------
build_go_figure <- function(ont, exclude_patterns, top_n = 10) {
  go_ont <- go_df %>%
    filter(ONTOLOGY == ont) %>%
    filter(!grepl(paste(exclude_patterns, collapse = "|"),
                  Description, ignore.case = TRUE)) %>%
    group_by(GeneSetTR) %>%
    slice_min(p.adjust, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      log10padj   = -log10(p.adjust),
      Description = str_wrap(Description, width = 50)
    )

  go_ont <- go_ont %>%
    group_by(GeneSetTR) %>%
    mutate(Description = forcats::fct_reorder(Description, log10padj)) %>%
    ungroup()

  p <- ggplot(go_ont, aes(x = log10padj, y = Description, fill = log10padj)) +
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

  n_terms <- length(unique(go_ont$Description))
  save_figure(p, paste0("Fig12_GO_", ont, "_enrichment_v3"), suffix = "",
              width = 8, height = max(5, n_terms * 0.33 + 2))
}

build_go_figure("BP", BP_EXCLUDE, top_n = 10)
build_go_figure("CC", CC_EXCLUDE, top_n = 10)
build_go_figure("MF", MF_EXCLUDE, top_n = 10)
