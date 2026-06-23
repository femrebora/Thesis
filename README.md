# Thesis

## Repository structure

```
Thesis/
├── README.md                       # This file
├── .gitignore                      # Data protection rules
├── .Rprofile                       # renv auto-load
├── renv.lock                       # Locked R environment
├── shared/                         # Cross-module utilities
│   ├── R/                          # Shared R functions
│   └── plotting/                   # Shared plotting themes
├── config/                         # Global configuration
├── docs/                           # Documentation
│   ├── methods_and_thresholds.md   # Analysis parameter table
│   └── figure_table_mapping.md     # Figure/table → script map
├── data/                           # (git-ignored — raw input data)
├── tests/                          # Validation tests
├── examples/                       # Example outputs
├── tcga_gbm/                       # Module 1: TCGA-GBM expression & survival
├── metabolomics/                   # Module 2: A172 metabolomics profiling
└── crispr_screen/                  # Module 3: CRISPR/Cas9 metabolic screen
```

## Modules

### 1. TCGA-GBM (`tcga_gbm/`)
12-gene expression and survival analysis in TCGA-GBM (TCGA-GBM cohort).
- **Entry point**: `Rscript tcga_gbm/scripts/run_tcga_analysis.R`
- **Key analyses**: RNA-seq expression, median-split Cox PH survival,
  continuous-expression Cox robustness, expression boxplots

### 2. Metabolomics (`metabolomics/`)
LC-MS metabolomics profiling of A172-S (sensitive) vs A172-R (resistant)
glioblastoma cell lines (104 metabolites, n=4 per group).
- **Entry point**: `Rscript metabolomics/scripts/run_metabolomics_analysis.R`
- **Key analyses**: Median/Log/Pareto normalization, Welch's t-test,
  fold-change analysis, volcano plot, pathway enrichment (SMPDB/KEGG ORA)

### 3. CRISPR/Cas9 Screen (`crispr_screen/`)
Genome-wide CRISPR/Cas9 metabolic gene knockout screen (Sabatini library,
~3,000 genes) in A172-R vs A172-S TRAIL-resistant GBM models.
- **Entry point**: `Rscript crispr_screen/scripts/run_crispr_analysis.R`
- **Key analyses**: MAGeCK RRA gene-level scoring, volcano plot, ranked LFC,
  heatmap, barplot, QC metrics, sgRNA-level profiles, cross-comparison,
  summary table, KEGG/Reactome/GO enrichment

## Quick start

```bash
# 1. Clone the repository
git clone <repo-url>
cd Thesis

# 2. Restore R environment
R -e 'renv::restore()'

# 3. Run analysis modules
Rscript tcga_gbm/scripts/run_tcga_analysis.R
Rscript metabolomics/scripts/run_metabolomics_analysis.R
Rscript crispr_screen/scripts/run_crispr_analysis.R
```

## Reproducibility

- **R version**: 4.4.2 (2024-10-31)
- **Environment**: Managed via [`renv`](https://rstudio.github.io/renv/)
- **Key packages**: TCGAbiolinks, MetaboAnalystR (v4.0.0), MAGeCK (v0.5.9.5),
  ggplot2, survival, survminer, EnhancedVolcano, mixOmics
- **Path anchoring**: All scripts use `here::here()` — run from repo root

## Data availability

Restricted data (TCGA patient-level downloads, raw CRISPR count matrices,
LC-MS raw files) are excluded from this repository per institutional policy.
Processed results, gene lists, and metadata templates are included where
they contain no protected health information.

See `docs/methods_and_thresholds.md` for full reproducibility parameters.

## Citation

See `CITATION.cff` (if available) or cite the thesis directly.

## License

All rights reserved. Contact the author for permissions.
