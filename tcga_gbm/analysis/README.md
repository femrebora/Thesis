# tcga_gbm/analysis/ — Phase 1: data acquisition, expression and survival

`01_expression_survival.R` is the analysis that produces every table the TCGA
figure scripts consume. It was previously missing from this repository —
`run_tcga_analysis.R` referred to it as living "in the source workspace" — so
Phase 1 could not be reproduced from a clone. It is now here.

```
01_expression_survival.R
   ├─ TCGAbiolinks GDCquery → GDCdata/ cache (STAR-counts + clinical)
   ├─ log2(TPM + 1) expression matrix for the 12-gene panel
   ├─ median-split Cox PH survival, per gene
   ├─ robustness: continuous Cox, multivariable Cox (+ age, sex)
   └─ writes ../results/01_expression_survival/tables/Table01–Table04
                    │
                    ▼
   ../scripts/make_tcga_figures.R  →  Şekil 4.13 – 4.15
```

## Running it

The script does `setwd(SCRIPT_DIR)` and then `source("config.R")`, and writes to
`../results/`. Run it from this directory so those relative paths resolve:

```bash
cd tcga_gbm/analysis
Rscript 01_expression_survival.R          # Phase 1 — downloads from GDC, slow
cd ../..
Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2   # Phase 2 — figures
```

Phase 1 needs network access and downloads several GB into `GDCdata/`
(git-ignored). Phase 2 needs only the committed Table01 / Table02 / Table04
CSVs for the Cox forest figures (Şekil 4.13–4.15), so **those thesis figures can
be regenerated without re-downloading TCGA**. `Table03_sample_qc_summary.csv` is
written by Phase 1 but is not committed.

## Two files named `utils.R`

| File | Used by | Contents |
|---|---|---|
| `tcga_gbm/analysis/utils.R` | `01_expression_survival.R` (Phase 1) | survival/plot helpers for the analysis |
| `tcga_gbm/R/utils.R` | `scripts/make_tcga_figures.R` (Phase 2) | Turkish figure theme, gene panel, `save_fig()` |

They are not interchangeable. Phase 1 sources its own by bare filename after
`setwd()`, which is why it must be run from this directory.

`config.R` is the single source of truth for the 12-gene panel, the thesis
Table 4.2 LFC values, and the analysis thresholds (`MIN_PATIENTS_KM = 20`,
`MIN_EVENTS_PER_GROUP = 3`).

## exploratory/

Analyses that were run but whose results are **not** in the submitted thesis.
Kept for provenance; none of them feed a thesis figure.

| Script | What it does |
|---|---|
| `02_ppi_network.R` | STRING v12 PPI network and hub-gene identification |
| `03_drug_interactions.py` | DGIdb v5 GraphQL drug–gene interaction mining |
| `04_immune_infiltration.R` | marker-signature immune deconvolution (MCP-counter-like) |
| `05_mutation_landscape.R` | MC3 consensus MAF mutation landscape |

`05_mutation_landscape.R` expects the MC3 MAF, which is ~200 MB and controlled
access — not distributed here.
