# FDR < 0.05 figure set

Every CRISPR figure rebuilt under a Benjamini–Hochberg **FDR < 0.05** acceptance
criterion instead of the uncorrected p < 0.05 used in the thesis figures.
Collected here so the corrected set can be reviewed as a unit.

Regenerate with:

```bash
Rscript crispr_screen/scripts/run_crispr_analysis_fdr05.R
```

The same figure scripts produce both sets; the criterion is switched via
`THESIS_SIG_METRIC` / `THESIS_SIG_THRESH` (see `crispr_screen/R/utils.R`).
`sessionInfo_fdr05.txt` records the R environment of the run that produced these.

## What each figure shows at FDR < 0.05

| File | Thesis ref | Result at FDR < 0.05 |
|---|---|---|
| `Fig01_volcano_fdr05` | Şekil 4.4 | **0 genes.** y-axis is −Log₁₀ FDR; the threshold line sits at 1.30 and the highest point in the screen reaches 0.37 |
| `Fig02_rank_lfc_fdr05` | Şekil 4.5 | **0 genes** highlighted; full LFC distribution retained in grey |
| `Fig03_heatmap_fdr05` | Şekil 4.6 | **0 genes** — null-result panel (a heatmap needs ≥1 row) |
| `Fig04_barplot_fdr05` | Şekil 4.7 | **0 genes** — three labelled null panels |
| `Fig05_qc_metrics_fdr05` | Şekil 4.1 | identical to baseline — QC does not depend on the threshold |
| `Fig06_sgrna_profiles_rebuilt_fdr05` | Şekil 4.3 | identical to baseline — raw sgRNA counts |
| `Fig08_cross_comparison_fdr05` | Şekil 4.8 | identical to baseline — LFC correlation uses all genes |
| `Fig09_summary_table_fdr05` | (not in thesis) | **0 genes** — null-result panel |
| `Fig10_kegg_enrichment_v3_fdr05` | Şekil 4.11 | 4 pathways (all genes); resistant-arm panel **empties** |
| `Fig11_reactome_enrichment_v3_fdr05` | Şekil 4.12 | 7 pathways (all genes), 1 (resistant arm: urea cycle) |
| `Fig12_GO_BP_enrichment_v3_fdr05` | Şekil 4.9 | top 10 shown; 47 / 16 terms available |
| `Fig12_GO_MF_enrichment_v3_fdr05` | Şekil 4.10 | top 10 shown; 25 / 67 terms available |
| `Fig12_GO_CC_enrichment_v3_fdr05` | (not in thesis) | 1 term per panel |

## Why the gene-level figures are empty

No gene in this screen reaches FDR < 0.05 — or even FDR < 0.25 — in any of the
three comparisons. Smallest observed FDR: 0.425 (Rpost/R0), 0.428 (Spost/Rpost),
0.500 (Spost/S0). This holds across every MAGeCK run of the dataset, so it is a
property of the screen, not of one parameter choice.

The scripts render this as an explicitly labelled null-result panel rather than a
blank canvas, so an empty figure can never be mistaken for a rendering failure.

Full threshold sensitivity analysis: `docs/methods_and_thresholds.md`.
