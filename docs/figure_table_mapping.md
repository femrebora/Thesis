# Figure & Table Mapping

Every final thesis figure and table mapped to its source script, input, and
output. This document is the reproduction guide for thesis-linked artefacts.

**Canonical thesis figure numbering** is the section
[Thesis figure numbering](#thesis-figure-numbering-canonical) below.
Other headings that mention Şekil numbers must match that section.

## Map format

| Field | Description |
|-------|-------------|
| Thesis ref | Figure/table number in the final thesis |
| Description | What the figure/table shows |
| Script | Source R script (relative to repo root) |
| Input | Required input data files |
| Output | Generated file(s) in module figures/ or results/ |
| Key function | Primary function generating the output |
| Reproduction command | Exact command to regenerate |

---

## Thesis figure numbering (canonical)

Repository documentation historically referred to the submitted thesis PDF as
`F_Emre_Bora_YL_Tez.pdf`. That PDF is **not distributed** in this repository.
The numbering below is the repository’s explicitly marked final thesis mapping.

The `FigNN` script names are internal; this is how they land in the submitted
thesis:

| Thesis ref | Script output | Description |
|---|---|---|
| Şekil 4.1 | `Fig05_qc_metrics` | Sample-level QC metrics |
| Şekil 4.2 | — (manual image) | NGS PCR gel photographs |
| Şekil 4.3 | `Fig06_sgrna_profiles_rebuilt` | sgRNA-level profiles, selected candidate genes |
| Şekil 4.4 | `Fig01_volcano` | Differential sgRNA representation (volcano) |
| Şekil 4.5 | `Fig02_rank_lfc` | Ranked log2 fold change |
| Şekil 4.6 | `Fig03_heatmap` | LFC heatmap of significant shared genes |
| Şekil 4.7 | `Fig04_barplot` | Significant genes by direction |
| Şekil 4.8 | `Fig08_cross_comparison` | Pearson correlation across comparisons |
| Şekil 4.9 | `Fig12_GO_BP_enrichment_v3` | GO biological process enrichment |
| Şekil 4.10 | `Fig12_GO_MF_enrichment_v3` | GO molecular function enrichment |
| Şekil 4.11 | `Fig10_kegg_enrichment_v3` | KEGG pathway enrichment |
| Şekil 4.12 | `Fig11_reactome_enrichment_v3` | Reactome pathway enrichment |
| Şekil 4.13 | `tcga_cox_orman` (HR panel) | TCGA-GBM candidate-gene hazard ratios |
| Şekil 4.14 | `tcga_cox_orman` | Median-split Cox forest plot |
| Şekil 4.15 | `tcga_surekli_cox_orman` | Continuous-expression Cox sensitivity |

**Not used in the submitted thesis:** `Fig09_summary_table`,
`Fig12_GO_CC_enrichment_v3`, the non-`_v3` enrichment figures, the entire
`metabolomics/` module, CRISPR `_fdr05` variants, and upstream Python figures.

Phase 2 also writes `tcga_ifade_cubuk.*` (expression summary bar plot fallback).
That file is generated for completeness; it is **not** assigned a Şekil number
in the canonical table above.

Every CRISPR figure additionally has a `_fdr05` variant produced by
`run_crispr_analysis_fdr05.R` (BH FDR < 0.05 acceptance criterion) — see
`docs/methods_and_thresholds.md` § FDR sensitivity analysis.

---

## TCGA-GBM Figures (Şekil 4.13–4.15)

| Thesis ref | Description | Script | Input | Output | Key function | Reproduction |
|------------|-------------|--------|-------|--------|-------------|-------------|
| Şekil 4.13–4.14 | Median-split Cox PH forest / HR panel | `tcga_gbm/scripts/make_tcga_figures.R` | `Table02_survival_summary.csv` | `tcga_gbm/figures/tcga_cox_orman.{pdf,png,tiff}` | `ggplot(geom_point + geom_errorbarh)` | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |
| Şekil 4.15 | Continuous Cox sensitivity forest plot | `tcga_gbm/scripts/make_tcga_figures.R` | `Table04_survival_robustness.csv` | `tcga_gbm/figures/tcga_surekli_cox_orman.{pdf,png,tiff}` | `ggplot(geom_point + geom_errorbarh)` | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |
| — (generated, no Şekil in canonical map) | Expression summary bar plot (fallback when GDC cache absent) | `tcga_gbm/scripts/make_tcga_figures.R` | `Table01_expression_gene_summary.csv` | `tcga_gbm/figures/tcga_ifade_cubuk.{pdf,png,tiff}` | `ggplot(bar)` | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |
| Survival summary | Turkish survival summary CSV | `tcga_gbm/scripts/make_tcga_figures.R` | `Table02_survival_summary.csv` | `tcga_gbm/figures/tcga_sagkalim_ozet_tablo.csv` | `write.csv()` with TR interpretation | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |

### TCGA prerequisite tables

These tables are produced by `01_expression_survival.R` (Phase 1) before Phase 2:

| Table | Source script | Rows | Tracked in git? | Description |
|-------|--------------|------|-----------------|-------------|
| `Table01_expression_gene_summary.csv` | `01_expression_survival.R` §5 | 12 genes | **Yes** | Mean ± SD expression per gene |
| `Table02_survival_summary.csv` | `01_expression_survival.R` §6 | 12 genes | **Yes** | Median-split Cox PH results |
| `Table03_sample_qc_summary.csv` | `01_expression_survival.R` §7 | N samples | **No** (generated locally by Phase 1) | Sample-level QC metrics |
| `Table04_survival_robustness.csv` | `01_expression_survival.R` §6' | 12 genes | **Yes** | Continuous Cox + multivariable Cox |

Phase 2 Cox figures do not require Table03.

---

## Metabolomics figures (PROVENANCE ONLY — not in submitted thesis)

> **Scope:** Metabolomics results were removed from the submitted thesis
> (unpublished data). The table below maps code → outputs for provenance. These
> are **not** final thesis figures.

| Ref | Description | Script | Input | Output | Reproduction |
|-----|-------------|--------|-------|--------|--------------|
| Provenance | S vs R volcano | `metabolomics/scripts/metabolomics_figures.R` | `metabolomics_clean.csv` | `metabolomics/figures/volcano.*`, `volcano_enhanced.*` | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Provenance | Top 9 metabolite boxplots | `metabolomics/scripts/metabolomics_figures.R` | `data_normalized_final.csv` | `metabolomics/figures/boxplots_top9.*` | same |
| Provenance | Heatmap | `metabolomics/scripts/metabolomics_figures.R` | `data_normalized_final.csv` | `metabolomics/figures/heatmap.*` | same |
| Provenance | Metabolite class pie | `metabolomics/scripts/metabolomics_figures.R` | `metabolite_classification_results.csv` | `metabolomics/figures/metabolite_class_pie.*` | same |
| Provenance | Sig by class | `metabolomics/scripts/metabolomics_figures.R` | `metabolite_classification_results.csv` | `metabolomics/figures/significant_by_class.*` | same |
| Provenance | Enrichment dot / category | `metabolomics/scripts/metabolomics_enrichment_pathway.R` | ORA inputs | `enrichment_*.*` | same |
| Exploratory | PCA / PLS-DA | `metabolomics/scripts/metabolomics_figures.R` | normalized data | `pca_*.*`, `plsda_*.*` | EXPLORATORY |
| Integration | CRISPR + metabolomics panels | `metabolomics/scripts/metabolomics_integration_figures.R` | integrated CSVs | `integration_*.*` | same |

See `metabolomics/upstream/README.md` for known runner limitations.

---

## CRISPR/Cas9 Screen Figures (Şekil 4.1–4.12)

Restricted inputs (`count_table.txt`, `*.gene_summary.tsv`) are required to
regenerate these; they are not in the public clone. Committed baseline figures
are the reference artefacts.

| Internal | Thesis ref | Description | Script | Output |
|----------|------------|-------------|--------|--------|
| Fig 01 | Şekil 4.4 | Volcano | `01_volcano.R` | `Fig01_volcano.{pdf,png}` |
| Fig 02 | Şekil 4.5 | Ranked LFC | `02_rank_lfc.R` | `Fig02_rank_lfc.{pdf,png}` |
| Fig 03 | Şekil 4.6 | Heatmap | `03_heatmap.R` | `Fig03_heatmap.{pdf,png}` |
| Fig 04 | Şekil 4.7 | Barplot | `04_barplot.R` | `Fig04_barplot.{pdf,png}` |
| Fig 05 | Şekil 4.1 | QC metrics | `05_qc_metrics.R` | `Fig05_qc_metrics.{pdf,png}` |
| Fig 06 | Şekil 4.3 | sgRNA profiles | `06_sgrna_profiles.R` | `Fig06_sgrna_profiles_rebuilt.{pdf,png}` |
| Fig 08 | Şekil 4.8 | Cross-comparison | `08_cross_comparison.R` | `Fig08_cross_comparison.{pdf,png}` |
| Fig 09 | — (not in thesis) | Summary table | `09_summary_table.R` | `Fig09_summary_table.*` + `results/top12_significant_genes.csv` |
| Fig 10 | Şekil 4.11 | KEGG (`_v3`) | `10_kegg_enrichment.R` | `Fig10_kegg_enrichment_v3.{pdf,png}` |
| Fig 11 | Şekil 4.12 | Reactome (`_v3`) | `11_reactome_enrichment.R` | `Fig11_reactome_enrichment_v3.{pdf,png}` |
| Fig 12 BP/MF | Şekil 4.9–4.10 | GO BP / MF (`_v3`) | `12_go_enrichment.R` | `Fig12_GO_{BP,MF}_enrichment_v3.{pdf,png}` |
| Fig 12 CC | — (not in thesis) | GO CC (`_v3`) | `12_go_enrichment.R` | `Fig12_GO_CC_enrichment_v3.{pdf,png}` |

Reproduction (authorised local inputs):
`Rscript crispr_screen/scripts/run_crispr_analysis.R`

---

## Table index

| Thesis ref | Description | Source table | Module | Location |
|------------|-------------|-------------|--------|----------|
| Tablo 4.2 | CRISPR LFC values (12 genes) | `GENE_LFC_RESISTANT` vector | crispr_screen → tcga_gbm | `tcga_gbm/R/utils.R` / `tcga_gbm/analysis/config.R` |
| — | Expression summary | `Table01_expression_gene_summary.csv` | tcga_gbm | `tcga_gbm/results/01_expression_survival/tables/` |
| — | Survival summary | `Table02_survival_summary.csv` | tcga_gbm | same |
| — | Sample QC | `Table03_sample_qc_summary.csv` | tcga_gbm | generated by Phase 1; **not committed** |
| — | Survival robustness | `Table04_survival_robustness.csv` | tcga_gbm | same tables dir |
| — | Turkish survival summary | `tcga_sagkalim_ozet_tablo.csv` | tcga_gbm | `tcga_gbm/figures/` (generated by Phase 2) |
| — | CRISPR top-12 / summary / enrichment CSVs | see `crispr_screen/results/` | crispr_screen | committed |
| — | Metabolomics processed tables | see `metabolomics/data/` | metabolomics | provenance only |
