# crispr_screen/ — CRISPR/Cas9 metabolic-gene screen

Computational module for the targeted **Sabatini metabolic-gene** CRISPR/Cas9
screen analysed in the Master’s thesis. This is **not** a genome-wide screen;
enrichment backgrounds follow the library gene universe (see
`load_library_universe()` in `R/utils.R` when used).

## What is final vs exploratory

| Artefact | Status |
|----------|--------|
| `scripts/01–06`, `08`, `10–12` (esp. `_v3` / `_rebuilt` outputs) | **Final thesis figures** (Şekil 4.1, 4.3–4.12) — see `docs/figure_table_mapping.md` |
| `scripts/09_summary_table.R`, GO CC panel | Generated; **not** used in the submitted thesis |
| `run_crispr_analysis_fdr05.R` / `figures/fdr05/` | **FDR sensitivity** documentation set |
| `05b_qc_palette_test.R`, `07_replicate_corr.R` | Exploratory / not final |
| `upstream/scripts/05–08_*.py` | Historical / sensitivity / first-generation Python figures |
| `upstream/nf-core/` | Provenance of the nf-core/crisprseq v2.3.0 R1-only pass |

Committed baseline figures under `figures/` are **immutable reference artefacts**.

## Experimental design of record (sample labels)

| Condition | Meaning |
|-----------|---------|
| **S0** | Sensitive (A172-S) — Day 0 baseline |
| **Spost** | Sensitive — tumour-derived after in-vivo growth |
| **R0** | Resistant (A172-R) — Day 0 baseline |
| **Rpost** | Resistant — tumour-derived after in-vivo growth |

There was **no TRAIL treatment in this screen**. Older upstream prose that says
“Post-TRAIL” is stale relative to this design of record.

Comparisons (MAGeCK): `Rpost_vs_R0`, `Spost_vs_Rpost`, `Spost_vs_S0`.

## Significance criterion of record

Submitted gene-level figures use **nominal MAGeCK RRA p < 0.05**.

- No gene reaches BH FDR < 0.05 (minimum observed FDR ≈ 0.425).
- Thesis Methods/Discussion text mentioning FDR < 0.25 does **not** match the
  executable figure workflow; see `docs/methods_and_thresholds.md`.
- This module preserves the criterion that generated the figures.

`SIG_THRESH` / `THESIS_SIG_*` in `R/utils.R` configure metric/suffix for the
FDR sensitivity runner. Several figure scripts still contain a local
`pval_thresh <- 0.05` (identical to the thesis default of 0.05). That behaviour
is documented rather than refactored without byte-identical regeneration proof.

## Public vs full local reproduction

**Public clone:** inspect code, parameters, committed gene lists, enrichment
CSVs, and baseline figures. Regenerating CRISPR figures will fail until
restricted inputs are present.

**Required restricted inputs** (place under `data/`; schemas in
`data/README.md`):

- `count_table.txt`
- `Rpost_vs_R0.gene_summary.tsv`
- `Spost_vs_Rpost.gene_summary.tsv`
- `Spost_vs_S0.gene_summary.tsv`

```bash
Rscript tools/check_environment.R
Rscript crispr_screen/scripts/run_crispr_analysis.R
# optional:
Rscript crispr_screen/scripts/run_crispr_analysis_fdr05.R
```

Upstream FASTQ → MAGeCK (authorised data only):

```bash
cd crispr_screen/upstream/scripts
bash run_all.sh --dry-run
bash run_all.sh
```

Protocol: [`upstream/REPRODUCIBLE_PIPELINE.md`](upstream/REPRODUCIBLE_PIPELINE.md) ·
[`upstream/README.md`](upstream/README.md)

## Directory map

```
crispr_screen/
├── README.md                 # this file
├── data/README.md            # restricted input schemas + gene_lists/
├── R/utils.R                 # labels, palettes, loaders, SIG_* config
├── scripts/                  # thesis figure pipeline + runners
├── results/                  # enrichment / summary CSVs, sessionInfo
├── figures/                  # baseline + historical + fdr05/
├── metadata/                 # sample_metadata_template.csv
└── upstream/                 # FASTQ → MAGeCK provenance
```
