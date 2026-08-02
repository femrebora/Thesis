# Methods & Thresholds

Reproducibility parameters for every analysis step across all three modules.
This table serves as the single source of truth for statistical thresholds,
normalization methods, and package versions.

## Table format

| Field | Description |
|-------|-------------|
| Module | tcga_gbm / metabolomics / crispr_screen |
| Step | Analysis step within the module |
| Method | Statistical method or algorithm |
| Package | R/Bioconductor package (with version) |
| Threshold | Cutoff, FDR, or significance level |
| Rationale | Why this threshold was chosen |
| Evidence source | Thesis section, protocol rule, or literature |
| Script | Source script implementing the step |
| Output | Generated file(s) |

---

## TCGA-GBM Module

| Module | Step | Method | Package | Threshold | Rationale | Evidence source | Script | Output |
|--------|------|--------|---------|-----------|-----------|-----------------|--------|--------|
| tcga_gbm | Data acquisition | GDCquery (RNA-seq STAR-counts) | TCGAbiolinks | TCGA-GBM project, HTSeq-Counts workflow | Standard GDC data access | Thesis §3.10 | `01_expression_survival.R` §1 | `GDCdata/` cache |
| tcga_gbm | Expression transform | log2(TPM + 1) | base R | — | Standard RNA-seq normalization; avoids log(0) | Thesis §3.10 | `01_expression_survival.R` §2 | Expression matrix |
| tcga_gbm | Sample type parsing | substr(barcode, 14, 15) | base R | 01 = Primary Tumor, 11 = Normal | TCGA barcode convention | GDC documentation | `utils.R` | sample_type column |
| tcga_gbm | Survival grouping | Median-split expression | survival (3.x) | Median gene expression | Standard approach; avoids data-driven cutpoint optimism | Thesis §3.10, Methods | `01_expression_survival.R` §6 | KM strata |
| tcga_gbm | Minimum patients | n ≥ 20 matched patients | — | MIN_PATIENTS_KM = 20 | Ensures adequate power for KM estimation | config.R | `01_expression_survival.R` §6 | Filtered cohort |
| tcga_gbm | Minimum events | n ≥ 3 events per group | — | MIN_EVENTS_PER_GROUP = 3 | Minimum for meaningful KM curve | config.R | `01_expression_survival.R` §6 | Filtered genes |
| tcga_gbm | Univariate Cox PH | coxph(Surv(OS.time, OS.event) ~ group) | survival | — | Median-split Cox, primary analysis | Thesis §3.10 | `01_expression_survival.R` §6 | HR, 95% CI, p |
| tcga_gbm | PH assumption check | cox.zph(cox_fit) | survival | p < 0.05 → PH violated | Schoenfeld residuals test; flagged in output | Thesis §3.10 | `01_expression_survival.R` §6 | PH_assumption_p |
| tcga_gbm | Continuous Cox (robustness) | coxph(Surv(OS.time, OS.event) ~ expr) | survival | — | Sensitivity analysis; avoids arbitrary dichotomy | Thesis §4.x (Şekil 4.17) | `01_expression_survival.R` §6' | cont_HR_per_SD |
| tcga_gbm | Multivariable Cox | coxph(Surv(OS.time, OS.event) ~ expr + age + sex) | survival | — | Adjusted for clinical covariates | Thesis §3.10 | `01_expression_survival.R` §6' | Adjusted HR |
| tcga_gbm | Optimal cutpoint (EXPLORATORY) | surv_cutpoint() | survminer | — | EXPLORATORY ONLY — not primary; optimism-biased | Thesis §3.10 | `01_expression_survival.R` §6' | cutpoint_logrank_p_EXPLORATORY |
| tcga_gbm | Multiple testing correction | p.adjust(method = "BH") | stats | FDR < 0.05 | Benjamini-Hochberg FDR across 12 genes | config.R | `01_expression_survival.R` §6' | adjusted_p |
| tcga_gbm | Figure output | ggsave PDF + PNG + TIFF | ggplot2 | 600 dpi PNG/TIFF, LZW compression | Publication requirement | config.R | `make_tcga_figures.R` | tcga_*.{pdf,png,tiff} |
| tcga_gbm | Forest plot x-axis | Log-scale with Turkish comma decimals | scales | breaks = 0.5, 0.7, 1.0, 1.3, 1.7 | Turkish locale convention (comma decimal) | make_tcga_figures_TR.R | `make_tcga_figures.R` | tcga_cox_orman.* |

## Metabolomics Module

| Module | Step | Method | Package | Threshold | Rationale | Evidence source | Script | Output |
|--------|------|--------|---------|-----------|-----------|-----------------|--------|--------|
| metabolomics | Missing value imputation | ReplaceMin (1/5 min positive) | MetaboAnalystR v4.0.0 | — | Standard for LC-MS metabolomics (left-censored missing) | MetaboAnalyst documentation | `metabolomics_figures.R` | Imputed matrix |
| metabolomics | Feature filtering | IQR filter, >25% missing | MetaboAnalystR v4.0.0 | IQR, 25% missing | Remove unreliable features | MetaboAnalyst documentation | `metabolomics_figures.R` | Filtered features |
| metabolomics | Normalization | MedianNorm + LogNorm + ParetoNorm | MetaboAnalystR v4.0.0 | — | Gold-standard for metabolomics; Pareto reduces heteroscedasticity | van den Berg et al. 2006 | `metabolomics_figures.R` | data_normalized_final.csv |
| metabolomics | Fold-change | FC.Anal(fc.thresh=2.0) | MetaboAnalystR v4.0.0 | |log2FC| > 1 | Standard 2-fold cutoff | TÜBİTAK analysis | `metabolomics_figures.R` | fold_change.csv |
| metabolomics | T-test | Welch's t-test (unequal var) | MetaboAnalystR v4.0.0 | p < 0.05 (nominal) | Small n (4 vs 4); Welch's robust to unequal variance | Thesis Methods | `metabolomics_figures.R` | t_test.csv |
| metabolomics | Multiple testing | FDR (Benjamini-Hochberg) | MetaboAnalystR v4.0.0 | FDR < 0.05 | Standard FDR correction | MetaboAnalyst documentation | `metabolomics_figures.R` | ttest_with_foldchange.csv |
| metabolomics | Volcano plot | log2FC vs -log10(p) | EnhancedVolcano | |log2FC| > 1, p < 0.05 | Standard volcano visualization | Thesis §4.x | `metabolomics_figures.R` | volcano.* |
| metabolomics | PCA (EXPLORATORY) | PCA.Anal() | MetaboAnalystR v4.0.0 | — | EXPLORATORY — excluded from thesis | Internal QC | `metabolomics_figures.R` | pca_score2d.* |
| metabolomics | PLS-DA (EXPLORATORY) | PLSR.Anal(reg=TRUE) | MetaboAnalystR v4.0.0 | — | EXPLORATORY — excluded from thesis | Internal QC | `metabolomics_figures.R` | plsda_*.* |
| metabolomics | Pathway enrichment (ORA) | Hypergeometric test | MetaboAnalystR v4.0.0 | FDR < 0.05 | Over-representation analysis against SMPDB/KEGG | TÜBİTAK reassessment | `metabolomics_enrichment_pathway.R` | enrichment_dot_plot.* |
| metabolomics | Metabolite classification | HMDB chemical class assignment | MetaboAnalystR v4.0.0 | — | Chemical taxonomy of significant metabolites | Internal | `metabolomics_figures.R` | metabolite_class_pie.* |

## CRISPR/Cas9 Screen Module

| Module | Step | Method | Package | Threshold | Rationale | Evidence source | Script | Output |
|--------|------|--------|---------|-----------|-----------|-----------------|--------|--------|
| crispr_screen | sgRNA counting | nf-core/crisprseq pipeline | Nextflow | — | Standardised CRISPR counting workflow | Thesis §3.6 | External (nf-core) | count_table.txt |
| crispr_screen | Gene-level scoring | MAGeCK RRA (v0.5.9.5) | MAGeCK | **Nominal p < 0.05 (UNCORRECTED)** | Robust Rank Aggregation. No gene reaches FDR < 0.05 (or even FDR < 0.25) in any comparison — smallest observed FDR = 0.425 — so the thesis figures report the uncorrected p-value. See §FDR sensitivity below. | Li et al. 2014 | `R/utils.R` `load_gene_summary()` | *.gene_summary.txt |
| crispr_screen | Direction classification | log2(LFC) sign concordance | utils.R | Concordant: same sign across all sig comparisons | Robust to single-comparison noise | utils.R | `01_volcano.R` | direction column |
| crispr_screen | Volcano plot | log2(Rpost/R0) vs -log10(p) | ggplot2 + ggrepel | Nominal p < 0.05, direction-concordant | y-axis is the UNCORRECTED p-value; no |LFC| cutoff is applied | Thesis Şekil 4.4 | `01_volcano.R` | Fig01_volcano.* |
| crispr_screen | Ranked LFC | Horizontal barplot, ranked by LFC | ggplot2 | — | Visualizes magnitude and direction | Thesis §4.x | `02_rank_lfc.R` | Fig02_rank_lfc.* |
| crispr_screen | Heatmap | ComplexHeatmap with NA preservation | ComplexHeatmap | Z-score scaling | NA cells preserved (not forced to 0) | Thesis §4.x | `03_heatmap.R` | Fig03_heatmap.* |
| crispr_screen | Barplot | Faceted by direction (Concordant classifier) | ggplot2 | — | Visualizes top significant genes | Thesis §4.x | `04_barplot.R` | Fig04_barplot.* |
| crispr_screen | QC metrics | Gini coefficient + Lorenz curve | ineq / base R | — | sgRNA-level evenness QC | Thesis §3.6-3.7 | `05_qc_metrics.R` | Fig05_qc_metrics.* |
| crispr_screen | sgRNA profiles | Per-gene sgRNA-level barplots | ggplot2 | Shared y-axis, IMP1 reference style | Individual sgRNA-level inspection | Thesis §4.x | `06_sgrna_profiles.R` | Fig06_sgrna_profiles.* |
| crispr_screen | Cross-comparison | Scatter plot Rpost vs Spost LFC | ggplot2 | — | Comparison across resistant/sensitive arms | Thesis §4.x | `08_cross_comparison.R` | Fig08_cross_comparison.* |
| crispr_screen | Summary table | Top-12 significant genes table | gridExtra / openxlsx | Top 12 by minimum nominal p across comparisons | Tabular summary for thesis | Thesis §4.x | `09_summary_table.R` | Fig09_summary_table.* |
| crispr_screen | KEGG enrichment | Over-representation analysis, BH-adjusted | clusterProfiler | **Top 10 by p.adjust, NO significance cutoff** | Panels may contain terms with p.adjust > 0.05; only 4/31 metabolic KEGG terms reach FDR < 0.05 | Thesis Şekil 4.11 | `10_kegg_enrichment.R` | Fig10_kegg_enrichment_v3.* |
| crispr_screen | Reactome enrichment | Over-representation analysis, BH-adjusted | ReactomePA | **Top 10 by p.adjust, NO significance cutoff** | 7/48 (all genes) and 1/33 (Rpost) metabolic terms reach FDR < 0.05 | Thesis Şekil 4.12 | `11_reactome_enrichment.R` | Fig11_reactome_enrichment_v3.* |
| crispr_screen | GO enrichment | Over-representation analysis, BH-adjusted | clusterProfiler | **Top 10 by p.adjust, NO significance cutoff**, BP+MF+CC | GO CC yields only 1 term per panel at FDR < 0.05 | Thesis Şekil 4.9-4.10 | `12_go_enrichment.R` | Fig12_GO_*_enrichment_v3.* |

## Reproducibility note: enrichment tables and KEGG drift

`crispr_screen/scripts/00_compute_enrichment.R` recomputes the six enrichment
tables from the two gene lists, with `pvalueCutoff = 0.1` and BH adjustment and
no `universe` argument. Re-running it was checked term by term against the
committed tables:

| Table | terms | IDs match | statistics match |
|---|---|---|---|
| `reactome_all_25.csv` | 95 | yes | **exact** |
| `reactome_Rpost_18.csv` | 73 | yes | **exact** |
| `go_all_25.csv` | 311 | yes | **exact** |
| `go_Rpost_18.csv` | 370 | yes | **exact** |
| `kegg_all_25.csv` | 86 | yes | small drift, see below |
| `kegg_Rpost_18.csv` | 73 | yes | small drift, see below |

**KEGG is a live database.** `enrichKEGG` fetches the current `hsa` annotation at
run time, and the human background has grown from **9,399 to 9,421** annotated
genes since the tables were first produced. For all 86 terms the gene membership
(`geneID`) and `GeneRatio` are unchanged — only the background denominator moved,
shifting p-values in the third to fourth decimal (e.g. hsa00983 p = 8.745e-4 →
8.686e-4). Term ranking and every downstream conclusion are unaffected.

The committed tables hold the values reported in the thesis and should not be
overwritten casually. To pin KEGG for an exact rerun, use a `KEGG.db`-style
snapshot or `clusterProfiler::download_KEGG()` cached to a fixed release.

## Reproducibility note: PDF vs PNG

`cairo_pdf` embeds a creation timestamp, and `save_as_docx` writes a zip with
per-run timestamps. Re-running any pipeline therefore always shows the `.pdf`
and `.docx` outputs as modified in `git status`, even when nothing about the
figure changed. **The `.png` outputs are byte-stable** and are the correct
artefact to diff when checking reproducibility:

```bash
Rscript crispr_screen/scripts/run_crispr_analysis.R
git status --short -- '*.png'    # empty == fully reproduced
```

## FDR sensitivity analysis (CRISPR module)

The thesis figures call a gene significant on its **uncorrected** MAGeCK RRA
p-value. To document what a conventional multiple-testing correction would do,
the entire figure set can be rebuilt under a Benjamini-Hochberg **FDR < 0.05**
acceptance criterion:

```bash
Rscript crispr_screen/scripts/run_crispr_analysis_fdr05.R
```

This sets `THESIS_SIG_METRIC=fdr`, `THESIS_SIG_THRESH=0.05` and
`THESIS_SIG_SUFFIX=_fdr05`; `R/utils.R` swaps the significance statistic and every
output is written with a `_fdr05` suffix, so no baseline figure is overwritten.
The identical figure scripts are used in both modes — there is no forked code path.

### Gene-level hit counts by threshold

| Comparison | genes tested | nominal p < 0.05 | FDR < 0.25 | FDR < 0.05 | smallest FDR |
|---|---|---|---|---|---|
| Rpost_vs_R0 (Dirençli: Tümör / Gün 0) | 178 | 18 | 0 | **0** | 0.425 |
| Spost_vs_Rpost (Duyarlı / Dirençli) | 164 | 15 | 0 | **0** | 0.428 |
| Spost_vs_S0 (Duyarlı: Tümör / Gün 0) | 13 | 2 | 0 | **0** | 0.500 |

No gene survives BH correction at any conventional threshold. Under `_fdr05` the
gene-level figures (Fig01–Fig04, Fig09) therefore render an explicitly labelled
null-result panel rather than a blank canvas. Fig05, Fig06 and Fig08 are
threshold-independent (QC, sgRNA counts, LFC correlation) and are unchanged.

### Enrichment terms surviving FDR < 0.05 (after the metabolic-term filter)

| Figure | Panel | terms shown at baseline | terms at FDR < 0.05 |
|---|---|---|---|
| Fig10 KEGG | Tüm genler | 10 | 4 |
| Fig10 KEGG | Dirençli: Tümör / Gün 0 | 10 | **0 (empty panel)** |
| Fig11 Reactome | Tüm genler | 10 | 7 |
| Fig11 Reactome | Dirençli: Tümör / Gün 0 | 10 | 1 |
| Fig12 GO BP | Tüm genler / Rpost | 10 / 10 | 10 / 10 (47 / 16 available) |
| Fig12 GO MF | Tüm genler / Rpost | 10 / 10 | 10 / 10 (25 / 67 available) |
| Fig12 GO CC | Tüm genler / Rpost | 1 / 10 | 1 / 1 |

Enrichment p-values are BH-adjusted by clusterProfiler/ReactomePA regardless of
mode; the `_fdr05` run adds the hard `p.adjust < 0.05` cutoff that the baseline
figures deliberately omit.

## Cross-module notes

### Palette conventions
- **CRISPR module**: Depleted (sgRNA loss) = blue (#2C7FB8), Enriched (sgRNA gain) = orange (#D95F02)
- **TCGA module**: Depleted (low expression) = vermillion (#D55E00), Enriched (high expression) = blue (#0072B2)
- **Rationale**: CRISPR and TCGA "Enriched/Depleted" have different biological meanings — CRISPR measures sgRNA representation change (enrichment = more sgRNA = gene loss phenotype), while TCGA measures mRNA abundance (enriched = higher expression). The inverse palettes reflect this. See thesis §4.x for full explanation.

### LFC value provenance
- CRISPR `GENE_LFC_RESISTANT` values (thesis Table 4.2) originate from the MAGeCK RRA Rpost_vs_R0 comparison
- These values are reused in TCGA visualizations as reference annotations
- See `tcga_gbm/R/utils.R` for the canonical LFC vector

### sgRNA library
- **Library**: Sabatini Human Metabolic Gene Knockout Library
- **File**: `library_sabatini_metko.tsv` (1.1 MB, ~3,000 genes)
- **Location**: Referenced from source workspace; not stored in this repo
- **Note**: Some genes have n=1 sgRNA, limiting statistical power for those targets

### Pipeline consistency
- All scripts use `here::here()` for path resolution — run from repository root
- All stochastic procedures use `set.seed(42)` for reproducibility
- All figure scripts write `sessionInfo()` output for full computational provenance
