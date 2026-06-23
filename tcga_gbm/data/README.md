# tcga_gbm/data/ — TCGA-GBM input data

## Directory structure

```
data/
├── README.md          # This file
└── (empty)            # All TCGA data is fetched at runtime via TCGAbiolinks
```

## Data sources

TCGA-GBM RNA-seq data is acquired programmatically at runtime via the
`TCGAbiolinks` Bioconductor package. No raw data files are stored in this
directory.

### Runtime downloads

| Data | Source | Method | Cache location |
|------|--------|--------|----------------|
| Gene expression (STAR-counts) | GDC Data Portal | `GDCquery()` | `../cache/GDCdata/` |
| Clinical data | GDC Data Portal | `GDCquery_clinic()` | `../cache/GDCdata/` |

### Cache

Downloaded GDC data is cached in `tcga_gbm/cache/GDCdata/` to avoid
re-downloading on every run. The cache is git-ignored.

## Precomputed tables

The survival analysis produces intermediate CSV tables stored in
`tcga_gbm/results/01_expression_survival/tables/`. These are the
inputs consumed by `make_tcga_figures.R`:

| Table | Description | Rows |
|-------|-------------|------|
| `Table01_expression_gene_summary.csv` | Per-gene mean ± SD expression | 12 genes |
| `Table02_survival_summary.csv` | Median-split Cox PH results | 12 genes |
| `Table03_sample_qc_summary.csv` | Sample-level QC metrics | N samples |
| `Table04_survival_robustness.csv` | Continuous Cox robustness | 12 genes |

## Gene panel

The 12-gene candidate panel is defined in `tcga_gbm/R/utils.R`:

- **Depleted** (CRISPR sgRNA depletion in resistant tumours): ASL, GBA3, ATP6V0A4, TPCN1
- **Enriched** (CRISPR sgRNA enrichment in resistant tumours): PLA2G4E, SLC1A1, RRM2, ARSD, GRIN1
- **Other** (not statistically significant in CRISPR): NAT2, SLC25A20, ALOX15

## .gitignore notes

- `cache/GDCdata/` is git-ignored (large TCGA downloads)
- `*.tsv`, `*.csv`, `*.xlsx` in `results/` are git-ignored except templates
