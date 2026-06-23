# Figure & Table Mapping

Every final thesis figure and table mapped to its source script, input, and
output. This document serves as the reproduction guide — to regenerate any
figure, find the script and run it from the repository root.

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

## TCGA-GBM Figures (Şekil 4.15–4.17)

| Thesis ref | Description | Script | Input | Output | Key function | Reproduction |
|------------|-------------|--------|-------|--------|-------------|-------------|
| Şekil 4.15 | 12-gene expression boxplot in TCGA-GBM | `tcga_gbm/scripts/make_tcga_figures.R` | `GDCdata/` cache (STAR-counts TSVs) or `Table01_expression_gene_summary.csv` | `tcga_gbm/figures/tcga_ifade_kutu.{pdf,png,tiff}` (or `tcga_ifade_cubuk.*` fallback) | `build_expr_from_cache()` → `ggplot(boxplot)` | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |
| Şekil 4.16 | Median-split Cox PH forest plot | `tcga_gbm/scripts/make_tcga_figures.R` | `Table02_survival_summary.csv` | `tcga_gbm/figures/tcga_cox_orman.{pdf,png,tiff}` | `ggplot(geom_point + geom_errorbarh)` with log10 x-axis | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |
| Şekil 4.17 | Continuous Cox sensitivity forest plot | `tcga_gbm/scripts/make_tcga_figures.R` | `Table04_survival_robustness.csv` | `tcga_gbm/figures/tcga_surekli_cox_orman.{pdf,png,tiff}` | `ggplot(geom_point + geom_errorbarh)` with FDR stars | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |
| Survival summary | Turkish survival summary CSV | `tcga_gbm/scripts/make_tcga_figures.R` | `Table02_survival_summary.csv` | `tcga_gbm/figures/tcga_sagkalim_ozet_tablo.csv` | `write.csv()` with TR interpretation | `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2` |

### TCGA prerequisite tables

These tables are produced by `01_expression_survival.R` (Phase 1) before Phase 2:

| Table | Source script | Rows | Description |
|-------|--------------|------|-------------|
| `Table01_expression_gene_summary.csv` | `01_expression_survival.R` §5 | 12 genes | Mean ± SD expression per gene |
| `Table02_survival_summary.csv` | `01_expression_survival.R` §6 | 12 genes | Median-split Cox PH results |
| `Table03_sample_qc_summary.csv` | `01_expression_survival.R` §7 | N samples | Sample-level QC metrics |
| `Table04_survival_robustness.csv` | `01_expression_survival.R` §6' | 12 genes | Continuous Cox + multivariable Cox |

---

## Metabolomics Figures

| Thesis ref | Description | Script | Input | Output | Key function | Reproduction |
|------------|-------------|--------|-------|--------|-------------|-------------|
| Volcano plot | S vs R volcano (log2FC vs -log10 p) | `metabolomics/scripts/metabolomics_figures.R` | `metabolomics_clean.csv` | `metabolomics/figures/volcano.{pdf,png,tiff}`, `volcano_enhanced.{pdf,png,tiff}` | EnhancedVolcano | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Boxplots (top 9) | Top 9 significant metabolites boxplots | `metabolomics/scripts/metabolomics_figures.R` | `data_normalized_final.csv` | `metabolomics/figures/boxplots_top9.{pdf,png,tiff}` | ggplot2 boxplot | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Heatmap | Hierarchical clustering heatmap | `metabolomics/scripts/metabolomics_figures.R` | `data_normalized_final.csv` | `metabolomics/figures/heatmap.{pdf,png,tiff}` | pheatmap | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Metabolite class pie | HMDB chemical class distribution | `metabolomics/scripts/metabolomics_figures.R` | `metabolite_classification_results.csv` | `metabolomics/figures/metabolite_class_pie.{pdf,png,tiff}` | ggplot2 pie | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Sig by class | Significant metabolites by chemical class | `metabolomics/scripts/metabolomics_figures.R` | `metabolite_classification_results.csv` | `metabolomics/figures/significant_by_class.{pdf,png,tiff}` | ggplot2 bar | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Enrichment dot plot | Pathway enrichment dot plot (ORA) | `metabolomics/scripts/metabolomics_enrichment_pathway.R` | `sig_metabolites_ora.csv` | `metabolomics/figures/enrichment_dot_plot.{pdf,png,tiff}` | ggplot2 dotplot | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Enrichment category bar | Enrichment by category | `metabolomics/scripts/metabolomics_enrichment_pathway.R` | `enrichment_category_results.csv` | `metabolomics/figures/enrichment_category_bar.{pdf,png,tiff}` | ggplot2 bar | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| PCA (EXPLORATORY) | PCA score plot (not in thesis) | `metabolomics/scripts/metabolomics_figures.R` | `data_normalized_final.csv` | `metabolomics/figures/pca_score2d.{pdf,png,tiff}`, `pca_scree.*` | MetaboAnalystR PCA | EXPLORATORY — excluded from thesis |
| PLS-DA (EXPLORATORY) | PLS-DA score + VIP (not in thesis) | `metabolomics/scripts/metabolomics_figures.R` | `data_normalized_final.csv` | `metabolomics/figures/plsda_*.*` | mixOmics PLS-DA | EXPLORATORY — excluded from thesis |

## Integration Figures (CRISPR + Metabolomics)

| Thesis ref | Description | Script | Input | Output | Key function | Reproduction |
|------------|-------------|--------|-------|--------|-------------|-------------|
| Dual volcano | Side-by-side: CRISPR genes + metabolites | `metabolomics/scripts/metabolomics_integration_figures.R` | `ttest_with_foldchange.csv`, CRISPR gene list | `metabolomics/figures/integration_dual_volcano.{pdf,png,tiff}` | ggplot2 + patchwork | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Pathway heatmap | Pathway coverage across omics layers | `metabolomics/scripts/metabolomics_integration_figures.R` | `integrated_pathway_analysis.csv` | `metabolomics/figures/integration_pathway_heatmap.{pdf,png,tiff}` | pheatmap | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Gene-metab network | Direct connections in shared pathways | `metabolomics/scripts/metabolomics_integration_figures.R` | `gene_metabolite_connections.csv` | `metabolomics/figures/integration_gene_metab_network.{pdf,png,tiff}` | ggplot2 network | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Rank comparison | Top hit rank comparison | `metabolomics/scripts/metabolomics_integration_figures.R` | `combined_metabolite_ranking.csv` | `metabolomics/figures/integration_rank_comparison.{pdf,png,tiff}` | ggplot2 rank plot | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |
| Pathway dotplot | Integrated pathway dot plot | `metabolomics/scripts/metabolomics_integration_figures.R` | `integrated_pathway_analysis.csv` | `metabolomics/figures/integration_pathway_dotplot.{pdf,png,tiff}` | ggplot2 dotplot | `Rscript metabolomics/scripts/run_metabolomics_analysis.R` |

## CRISPR/Cas9 Screen Figures (Şekil 4.1–4.12)

| Thesis ref | Description | Script | Input | Output | Key function | Reproduction |
|------------|-------------|--------|-------|--------|-------------|-------------|
| Fig 01 | Volcano plot (Diferansiyel sgRNA temsili) | `crispr_screen/scripts/01_volcano.R` | MAGeCK gene_summary.txt | `crispr_screen/figures/Fig01_volcano.{pdf,png}` | EnhancedVolcano | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 02 | Ranked LFC (Sıralı log2 kat değişimi) | `crispr_screen/scripts/02_rank_lfc.R` | MAGeCK gene_summary.txt + concordant classifier | `crispr_screen/figures/Fig02_rank_lfc.{pdf,png}` | ggplot2 bar | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 03 | Heatmap (Gen temsili ısı haritası) | `crispr_screen/scripts/03_heatmap.R` | MAGeCK gene_summary.txt (all comparisons) | `crispr_screen/figures/Fig03_heatmap.{pdf,png}` | ComplexHeatmap | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 04 | Barplot (Öncelikli genler) | `crispr_screen/scripts/04_barplot.R` | MAGeCK gene_summary.txt + concordant classifier | `crispr_screen/figures/Fig04_barplot.{pdf,png}` | ggplot2 bar (faceted) | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 05 | QC metrics (Kalite kontrol ölçütleri) | `crispr_screen/scripts/05_qc_metrics.R` | count_table.txt | `crispr_screen/figures/Fig05_qc_metrics.{pdf,png}` | Gini + Lorenz + ggplot2 | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 06 | sgRNA profiles (sgRNA düzeyinde profiller) | `crispr_screen/scripts/06_sgrna_profiles.R` | count_table.txt | `crispr_screen/figures/Fig06_sgrna_profiles.{pdf,png}` | ggplot2 bar (shared y) | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 08 | Cross-comparison (Çapraz karşılaştırma) | `crispr_screen/scripts/08_cross_comparison.R` | MAGeCK gene_summary.txt | `crispr_screen/figures/Fig08_cross_comparison.{pdf,png}` | ggplot2 scatter | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 09 | Summary table (Özet tablo) | `crispr_screen/scripts/09_summary_table.R` | MAGeCK gene_summary.txt | `crispr_screen/figures/Fig09_summary_table.{pdf,png}` + `results/top12_significant_genes.csv` | gridExtra table | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 10 | KEGG enrichment | `crispr_screen/scripts/10_kegg_enrichment.R` | gene lists + org.Hs.eg.db | `crispr_screen/figures/Fig10_kegg_enrichment.{pdf,png}` + `results/kegg_*.csv` | clusterProfiler enrichKEGG | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 11 | Reactome enrichment | `crispr_screen/scripts/11_reactome_enrichment.R` | gene lists + org.Hs.eg.db | `crispr_screen/figures/Fig11_reactome_enrichment.{pdf,png}` + `results/reactome_*.csv` | ReactomePA enrichPathway | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |
| Fig 12 | GO enrichment (BP + MF + CC) | `crispr_screen/scripts/12_go_enrichment.R` | gene lists + org.Hs.eg.db | `crispr_screen/figures/Fig12_GO_*_enrichment.{pdf,png}` + `results/go_*.csv` | clusterProfiler enrichGO | `Rscript crispr_screen/scripts/run_crispr_analysis.R` |

## Table index

| Thesis ref | Description | Source table | Module | Location |
|------------|-------------|-------------|--------|----------|
| Tablo 4.2 | CRISPR LFC values (12 genes) | `GENE_LFC_RESISTANT` vector | crispr_screen → tcga_gbm | `tcga_gbm/R/utils.R` |
| — | Expression summary | `Table01_expression_gene_summary.csv` | tcga_gbm | `tcga_gbm/results/01_expression_survival/tables/` |
| — | Survival summary | `Table02_survival_summary.csv` | tcga_gbm | `tcga_gbm/results/01_expression_survival/tables/` |
| — | Sample QC | `Table03_sample_qc_summary.csv` | tcga_gbm | `tcga_gbm/results/01_expression_survival/tables/` |
| — | Survival robustness | `Table04_survival_robustness.csv` | tcga_gbm | `tcga_gbm/results/01_expression_survival/tables/` |
| — | Turkish survival summary | `tcga_sagkalim_ozet_tablo.csv` | tcga_gbm | `tcga_gbm/figures/` |
| — | CRISPR top-12 genes | `top12_significant_genes.csv` | crispr_screen | `crispr_screen/results/` |
| — | CRISPR summary | `summary_significant_genes.csv` | crispr_screen | `crispr_screen/results/` |
| — | KEGG enrichment (all) | `kegg_all_25.csv` | crispr_screen | `crispr_screen/results/` |
| — | KEGG enrichment (Rpost) | `kegg_Rpost_18.csv` | crispr_screen | `crispr_screen/results/` |
| — | Reactome enrichment (all) | `reactome_all_25.csv` | crispr_screen | `crispr_screen/results/` |
| — | Reactome enrichment (Rpost) | `reactome_Rpost_18.csv` | crispr_screen | `crispr_screen/results/` |
| — | GO enrichment (all) | `go_all_25.csv` | crispr_screen | `crispr_screen/results/` |
| — | GO enrichment (Rpost) | `go_Rpost_18.csv` | crispr_screen | `crispr_screen/results/` |
| — | T-test + fold-change | `ttest_with_foldchange.csv` | metabolomics | `metabolomics/data/` |
| — | Metabolite classification | `metabolite_classification_results.csv` | metabolomics | `metabolomics/data/` |
| — | Pathway integration | `integrated_pathway_analysis.csv` | metabolomics | `metabolomics/data/` |
