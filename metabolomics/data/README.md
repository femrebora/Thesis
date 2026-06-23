# metabolomics/data/ — A172-S vs A172-R LC-MS metabolomics

## Study design

| Group | Cell line | TRAIL sensitivity | Replicates |
|-------|-----------|-------------------|------------|
| A172-S | A172      | Sensitive         | n = 4      |
| A172-R | A172      | Resistant         | n = 4      |

## Input data

| File | Description | Rows |
|------|-------------|------|
| `metabolomics_clean.csv` | Cleaned peak intensity table | 8 samples × 104 metabolites |
| `metabolomics_conc_format.csv` | Concentration-format data (colu layout) | 8 rows × 106 cols |
| `metabolomics_tidy.csv` | Tidy/long-format data | — |

## Processed data (script outputs)

| File | Description | Source script |
|------|-------------|---------------|
| `data_normalized_final.csv` | MedianNorm + LogNorm + ParetoNorm data | `metabolomics_figures.R` |
| `ttest_with_foldchange.csv` | T-test results with fold-change | `metabolomics_figures.R` |
| `fold_change.csv` | Log2 fold-change (S vs R) | `metabolomics_figures.R` |
| `sig_metabolites_ora.csv` | Significant metabolites for ORA input | `metabolomics_figures.R` |
| `metabolite_classification_results.csv` | HMDB chemical class assignments | `metabolomics_figures.R` |
| `metabolite_hmdb_mapping.csv` | Metabolite → HMDB ID mapping | `metabolomics_enrichment_pathway.R` |
| `enrichment_category_results.csv` | Category-level enrichment results | `metabolomics_enrichment_pathway.R` |
| `plsda_vip_scores.csv` | PLS-DA VIP scores (exploratory only) | `metabolomics_figures.R` |
| `combined_metabolite_ranking.csv` | Multi-metric metabolite ranking | `metabolomics_integration_figures.R` |
| `integrated_pathway_analysis.csv` | CRISPR + metabolomics integrated pathway table | `metabolomics_integration_figures.R` |
| `gene_metabolite_connections.csv` | Direct gene-metabolite links | `metabolomics_integration_figures.R` |

## Metabolite identification

Metabolites were identified via:
1. Accurate mass matching against HMDB, METLIN, and LipidMaps
2. MS/MS spectral matching where available
3. Retention time matching against in-house standards

⚠️ **Note**: 2,6-Dichloro-4-nitroaniline (top hit) may be an environmental
contaminant — flagged for verification in the thesis.

## .gitignore notes

- `*.csv`, `*.tsv`, `*.xlsx` in `data/` are git-ignored except README.md
- Raw LC-MS files (.raw, .mzML, .mzXML) are excluded
