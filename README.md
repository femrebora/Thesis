# Thesis

> **Scope note.** The metabolomics module is retained here for provenance but its
> results were removed from the submitted thesis (unpublished data). The thesis
> figure set is produced by the `crispr_screen` and `tcga_gbm` modules only.

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
│   ├── analysis/                   #   Phase 1 — GDC download, expression, Cox survival
│   │   └── exploratory/            #   PPI, drug, immune, mutation (not in thesis)
│   ├── scripts/                    #   Phase 2 — figures from precomputed tables
│   └── results/                    #   Table01–Table04 (committed)
├── metabolomics/                   # Module 2: A172 metabolomics profiling
│   ├── upstream/                   #   MetaboAnalystR session records (provenance)
│   └── scripts/                    #   figure code
└── crispr_screen/                  # Module 3: CRISPR/Cas9 metabolic screen
    ├── upstream/                   #   FASTQ → MAGeCK count → MAGeCK test
    │   ├── nf-core/                #   nf-core/crisprseq v2.3.0 run
    │   └── scripts/                #   01–08 + run_all.sh (pipeline of record)
    ├── scripts/                    #   00 enrichment + 01–12 thesis figures
    └── results/                    #   enrichment tables, gene summaries
```

## Analysis chain, end to end

```
CRISPR   FASTQ ─▶ upstream/scripts/run_all.sh ─▶ count_table.txt
                                              ─▶ *.gene_summary.tsv
                     ─▶ scripts/00_compute_enrichment.R ─▶ results/{kegg,reactome,go}_*.csv
                     ─▶ scripts/run_crispr_analysis.R   ─▶ Şekil 4.1, 4.3–4.12

TCGA     GDC ─▶ analysis/01_expression_survival.R ─▶ results/.../Table01–04
                     ─▶ scripts/run_tcga_analysis.R phase2 ─▶ Şekil 4.13–4.15
```

Each stage is independently runnable from the artefacts the previous one commits,
so the thesis figures can be rebuilt without re-running the restricted upstream
steps or re-downloading TCGA.

## Modules

### 1. TCGA-GBM (`tcga_gbm/`)
12-gene expression and survival analysis in the TCGA-GBM cohort.
- **Phase 1 (data + statistics)**: `cd tcga_gbm/analysis && Rscript 01_expression_survival.R`
  — downloads from GDC, writes Table01–Table04. Needs network; slow.
- **Phase 2 (figures)**: `Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2`
  — reads the committed tables, no download needed.
- **Key analyses**: RNA-seq expression, median-split Cox PH survival,
  continuous-expression Cox robustness, expression boxplots

### 2. Metabolomics (`metabolomics/`)
LC-MS metabolomics profiling of A172-S (sensitive) vs A172-R (resistant)
glioblastoma cell lines (104 metabolites, n=4 per group).
- **Entry point**: `Rscript metabolomics/scripts/run_metabolomics_analysis.R`
  (known issue: does not self-terminate — see `metabolomics/upstream/README.md`)
- **Key analyses**: Median/Log/Pareto normalization, Welch's t-test,
  fold-change analysis, volcano plot, pathway enrichment (SMPDB/KEGG ORA)

### 3. CRISPR/Cas9 Screen (`crispr_screen/`)
Targeted CRISPR/Cas9 **metabolic-gene** knockout screen (Sabatini metabolic
library) in A172-R vs A172-S TRAIL-resistant GBM models. This is a focused
sub-library screen, not a genome-wide screen; the enrichment background is the
set of genes the library targets (see `load_library_universe()` in
`crispr_screen/R/utils.R`).
- **Upstream (restricted inputs)**: `cd crispr_screen/upstream/scripts && bash run_all.sh`
  — FASTQ → MAGeCK count → MAGeCK test. See `crispr_screen/upstream/README.md`.
- **Enrichment tables**: `Rscript crispr_screen/scripts/00_compute_enrichment.R`
  — needs network (KEGG REST). Only needed if regenerating the tables.
- **Figures**: `Rscript crispr_screen/scripts/run_crispr_analysis.R`
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

# 3. Regenerate the thesis figures from the committed intermediate artefacts.
#    Neither step needs the restricted raw data or a TCGA download.
Rscript crispr_screen/scripts/run_crispr_analysis.R    # Şekil 4.1, 4.3–4.12
Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2    # Şekil 4.13–4.15

# 4. (optional) Rebuild the CRISPR figures under a BH FDR < 0.05 criterion.
#    Outputs carry a "_fdr05" suffix; baseline figures are never overwritten.
Rscript crispr_screen/scripts/run_crispr_analysis_fdr05.R
```

## Significance criterion

The threshold is defined once, as `SIG_THRESH` in `crispr_screen/R/utils.R`.
To change it, re-derive the gene lists and enrichment tables as well:

```bash
export THESIS_SIG_THRESH=0.01
Rscript crispr_screen/scripts/00_derive_gene_lists.R
Rscript crispr_screen/scripts/00_compute_enrichment.R
Rscript crispr_screen/scripts/run_crispr_analysis.R
```

`run_crispr_analysis.R` refuses to render if the gene lists on disk do not match
the threshold in force, so a partially-updated figure set cannot be produced
silently. Verify at any time with
`Rscript crispr_screen/scripts/00_derive_gene_lists.R --check`.

Gene-level CRISPR significance in the thesis figures is the **uncorrected**
MAGeCK RRA p-value (p < 0.05). No gene in this screen reaches FDR < 0.05 in any
comparison (smallest observed FDR = 0.425). `docs/methods_and_thresholds.md`
documents the full threshold sensitivity analysis and how to reproduce the
FDR-corrected figure set.

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
