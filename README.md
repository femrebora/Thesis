# Master’s Thesis — Computational Analysis Repository

**Author:** Furkan Emre Bora
**Repository:** [https://github.com/femrebora/Thesis](https://github.com/femrebora/Thesis)

This repository contains the computational analysis, reproducibility materials, processed outputs, and figure-generation code associated with the author’s Master’s thesis.

It is a **research compendium / thesis reproducibility archive**, not a general-purpose software project. The scientific analysis that produced the submitted thesis figures is treated as a **frozen record**: parameters, gene lists, statistical thresholds, sample definitions, and committed baseline figures are preserved as used for the thesis.

> The thesis PDF itself is **not distributed** in this repository. Documentation refers to the submitted thesis by the working filename `F_Emre_Bora_YL_Tez.pdf` where that filename appears in historical notes; the file is not present in this clone.

---

## Thesis relationship

This GitHub repository is the computational companion to the submitted Master’s thesis. It documents:

- which code generated which thesis figure;
- which parameters and software versions were used;
- which inputs are public versus restricted;
- what was exploratory or excluded from the final thesis.

It does **not** re-analyse the screen under newer methodological preferences.

---

## Research overview

The thesis work centres on metabolic aspects of TRAIL resistance in glioblastoma (GBM) models, with three computational modules:

| Module | Role in the submitted thesis |
|--------|------------------------------|
| **CRISPR/Cas9 metabolic screen** (`crispr_screen/`) | Final thesis figures (Şekil 4.1, 4.3–4.12) |
| **TCGA-GBM** (`tcga_gbm/`) | Final thesis figures (Şekil 4.13–4.15) |
| **Metabolomics** (`metabolomics/`) | **Provenance only — not included in the submitted thesis results** |

---

## What this repository contains

- Executable R (and upstream Python) workflows used for the thesis analyses
- Documented analysis parameters and thresholds (`docs/methods_and_thresholds.md`)
- Figure/table → script mapping (`docs/figure_table_mapping.md`)
- Committed processed, non-sensitive tables and gene lists
- Committed baseline thesis figures (reference artefacts)
- Historical environment manifest (`ENVIRONMENT.txt`)
- Upstream CRISPR FASTQ → MAGeCK provenance
- Metabolomics session/figure code retained for research-process documentation

---

## Final thesis analysis scope

**Included in the submitted thesis figure set**

- CRISPR R figure pipeline (`crispr_screen/scripts/`), criterion of record: nominal MAGeCK RRA **p < 0.05**
- TCGA-GBM Phase 1 tables + Phase 2 Cox forest figures

**Retained but not part of the submitted thesis results**

- Entire `metabolomics/` module (unpublished / excluded from the final thesis)
- CRISPR `_fdr05` sensitivity figure set
- Upstream Python first-generation figures and exploratory CRISPR scripts (`05b`, `07_replicate_corr`, improved/scaffold analyses)
- TCGA `analysis/exploratory/` (PPI, drug, immune, mutation)

---

## Repository structure

```
Thesis/
├── README.md
├── ENVIRONMENT.txt                 # historical R/package versions (thesis run)
├── .Rprofile                       # loads renv only if renv.lock exists
├── .gitignore
├── tools/
│   └── check_environment.R         # validate installed versions vs ENVIRONMENT.txt
├── docs/
│   ├── methods_and_thresholds.md   # analysis parameters (source of truth)
│   └── figure_table_mapping.md     # thesis Şekil/Tablo → script map
├── shared/
│   └── R/plotting_themes.R         # shared colour helpers (optional utility)
├── crispr_screen/
│   ├── README.md                   # module overview + public vs local reproduction
│   ├── data/                       # restricted inputs (git-ignored) + gene_lists/
│   ├── R/utils.R                   # shared CRISPR helpers / labels / palettes
│   ├── scripts/                    # 00_* + 01–12 thesis figure scripts + runners
│   ├── results/                    # enrichment tables, summaries, sessionInfo
│   ├── figures/                    # baseline + _fdr05 / historical variants
│   ├── metadata/                   # sample metadata template
│   └── upstream/                   # FASTQ → MAGeCK pipeline of record
├── tcga_gbm/
│   ├── analysis/                   # Phase 1 (GDC + survival) + exploratory/
│   ├── scripts/                    # Phase 2 figures + module runner
│   ├── results/.../tables/         # Table01, Table02, Table04 (committed)
│   └── figures/                    # Şekil 4.13–4.15 outputs
└── metabolomics/                   # PROVENANCE / unpublished / not in thesis
    ├── upstream/                   # MetaboAnalystR session transcripts
    ├── scripts/                    # figure regeneration (known limitations)
    ├── data/                       # processed metabolite tables
    └── figures/
```

There is **no** top-level `config/`, `data/`, `tests/`, or `examples/` directory, and no `shared/plotting/` directory (plotting helpers live in `shared/R/`).

---

## Reproducibility levels

### 1. Public repository (fresh clone)

A public clone supports **inspection and partial reproduction**:

- exact code and documented parameters;
- environment version manifest;
- committed gene lists and processed non-sensitive CRISPR enrichment/summary tables;
- committed thesis figures (as reference artefacts);
- TCGA Phase 2 figures from committed summary tables (no GDC download required for the Cox forest panels);
- methodology documentation and figure provenance.

A fresh public clone **cannot** regenerate CRISPR thesis figures: the required count matrix and MAGeCK gene summaries are intentionally not distributed.

### 2. Full local / authorised reproduction

With restricted CRISPR inputs placed under `crispr_screen/data/` (schemas in `crispr_screen/data/README.md`), the author or an authorised researcher can regenerate CRISPR figures with the recorded environment.

TCGA Phase 1 additionally requires network access to GDC via `TCGAbiolinks`.

---

## Quick start

### Inspection / public clone

```bash
git clone https://github.com/femrebora/Thesis.git
cd Thesis

# Compare your R/packages to the historical thesis environment
Rscript tools/check_environment.R

# TCGA Phase 2 — regenerates Cox forest figures from committed tables
# (no private CRISPR inputs; no GDC download needed for those panels)
Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2
```

CRISPR figure regeneration on a fresh clone is expected to **stop with a clear error** listing the missing restricted files.

### Full authorised / local reproduction (restricted CRISPR inputs)

Place these files under `crispr_screen/data/` (see `crispr_screen/data/README.md`):

- `count_table.txt`
- `Rpost_vs_R0.gene_summary.tsv`
- `Spost_vs_Rpost.gene_summary.tsv`
- `Spost_vs_S0.gene_summary.tsv`

Then, ideally under R **4.4.2** with package versions matching `ENVIRONMENT.txt`:

```bash
Rscript tools/check_environment.R
Rscript crispr_screen/scripts/run_crispr_analysis.R          # Şekil 4.1, 4.3–4.12
Rscript tcga_gbm/scripts/run_tcga_analysis.R phase2          # Şekil 4.13–4.15

# Optional sensitivity set (writes *_fdr05 outputs; does not overwrite baseline)
Rscript crispr_screen/scripts/run_crispr_analysis_fdr05.R
```

**Network notes**

- TCGA Phase 1 (`tcga_gbm/analysis/01_expression_survival.R`) downloads from GDC.
- CRISPR enrichment recomputation (`00_compute_enrichment.R`) contacts live KEGG/Reactome resources; committed enrichment CSVs are the thesis values of record.
- Metabolomics runner may hang on network calls; see `metabolomics/upstream/README.md`.

Compare regenerated CRISPR figures via **PNG** digests (PDFs embed timestamps):

```bash
git status --short -- '*.png'
```

---

## CRISPR/Cas9 workflow

Targeted **metabolic-gene** knockout screen (Sabatini Human Metabolic Gene Knockout Library), not a genome-wide CRISPR screen.

| Layer | Location | Role |
|-------|----------|------|
| Upstream pipeline of record | `crispr_screen/upstream/scripts/run_all.sh` | FASTQ → MAGeCK count → MAGeCK test |
| nf-core provenance | `crispr_screen/upstream/nf-core/` | crisprseq v2.3.0 (R1-only first pass; superseded for counts) |
| MAGeCK intermediates | `crispr_screen/data/*.gene_summary.tsv`, `count_table.txt` | **Restricted — not in public repo** |
| Final R figure pipeline | `crispr_screen/scripts/run_crispr_analysis.R` | Submitted thesis figures |
| FDR sensitivity | `run_crispr_analysis_fdr05.R` + `figures/fdr05/` | Exploratory / documentation |
| Exploratory scripts | `05b_qc_palette_test.R`, `07_replicate_corr.R`, upstream `06–08` | Not final thesis figures |

**Criterion of record for submitted gene-level figures:** nominal MAGeCK RRA **p < 0.05**.
No gene reaches BH FDR < 0.05 (minimum observed FDR ≈ 0.425). Parts of the written thesis mention FDR < 0.25; applying that criterion would not reproduce the reported hit counts. This repository preserves the computational criterion that generated the figures and records the wording discrepancy in `docs/methods_and_thresholds.md`.

Post samples (**Spost**, **Rpost**) are tumour-derived after in-vivo growth; there was **no TRAIL treatment in this screen**. Resistance is a cell-line property (A172-S vs A172-R).

Module docs: [`crispr_screen/README.md`](crispr_screen/README.md) · [`crispr_screen/upstream/REPRODUCIBLE_PIPELINE.md`](crispr_screen/upstream/REPRODUCIBLE_PIPELINE.md)

---

## TCGA-GBM workflow

Public/controlled-access data originate from **TCGA / GDC** (project **TCGA-GBM**) and are subject to GDC access and use terms. They are not author-owned experimental data.

| Phase | Script | Needs | Output |
|-------|--------|-------|--------|
| **Phase 1** | `tcga_gbm/analysis/01_expression_survival.R` (also mirrored under `scripts/`) | Network + GDC download | Table01–Table04 under `results/.../tables/` |
| **Phase 2** | `tcga_gbm/scripts/make_tcga_figures.R` via `run_tcga_analysis.R phase2` | Committed tables | Şekil 4.13–4.15 |

Frozen scientific configuration includes the **12-gene panel**, thesis **Tablo 4.2** LFC values, median-split Cox, continuous Cox, multivariable models, and FDR settings in `tcga_gbm/analysis/config.R` / `tcga_gbm/R/utils.R`.

**Note:** `Table03_sample_qc_summary.csv` is produced by Phase 1 but is **not** committed in this repository. Table01, Table02, and Table04 are committed.

See [`tcga_gbm/analysis/README.md`](tcga_gbm/analysis/README.md).

---

## Metabolomics provenance note

**PROVENANCE / UNPUBLISHED / NOT INCLUDED IN THE FINAL SUBMITTED THESIS RESULTS**

The metabolomics module documents broader research process (LC-MS profiling of A172-S vs A172-R). Do not treat its plots as final thesis figures. Session files under `metabolomics/upstream/` are interactive MetaboAnalystR transcripts, not clean standalone scripts. Known runner/network limitations are documented in `metabolomics/upstream/README.md` and do not affect reproduction of the submitted thesis.

---

## Environment and software versions

`ENVIRONMENT.txt` is the **historical environment manifest** for the thesis figure run (R 4.4.2, captured 2026-08-02). Do not silently replace those versions with whatever is installed today.

```bash
Rscript tools/check_environment.R
```

There is **no** `renv.lock` in this repository. `.Rprofile` loads renv only if a lockfile is later added; do not run `renv::init()` expecting to recreate the thesis environment from today’s packages.

Key tools used in the broader workflow (see module docs for detail): MAGeCK 0.5.9.5, nf-core/crisprseq v2.3.0, TCGAbiolinks, clusterProfiler / ReactomePA, MetaboAnalystR v4.0.0 (metabolomics provenance).

---

## Analysis parameters / frozen-analysis statement

Detailed thresholds, comparisons, and package notes live in [`docs/methods_and_thresholds.md`](docs/methods_and_thresholds.md).

This archive preserves what was actually run for the submitted figures. Where thesis wording differs from the executable workflow, the discrepancy is documented rather than retroactively altering the analysis.

---

## Figure and table provenance

Complete mapping: [`docs/figure_table_mapping.md`](docs/figure_table_mapping.md).

Committed baseline PNGs/PDFs are **immutable reference artefacts**. Prefer documenting historical/`_v3`/`_rebuilt`/`_fdr05` variants over deleting them.

---

## Data availability and data-source statement

| Source | What | In this repo |
|--------|------|--------------|
| Author / laboratory experimental CRISPR sequencing | FASTQ, count matrix, MAGeCK gene summaries | **Not redistributed** publicly (git-ignored); schemas documented |
| Author / laboratory metabolomics | Raw LC-MS; processed metabolite tables | Raw not included; processed CSVs retained for provenance |
| TCGA / GDC | TCGA-GBM RNA-seq + clinical | Downloaded at runtime; large caches git-ignored; summary tables partly committed |
| Processed non-sensitive outputs | Gene lists, enrichment CSVs, TCGA Table01/02/04, figures | Distributed for inspection / partial reproduction |

Do not assume a single ownership model for all files. Third-party TCGA/GDC data remain subject to their access and use terms.

---

## Known limitations / reproducibility notes

- Public clone cannot regenerate CRISPR thesis figures without restricted inputs.
- `SIG_THRESH` in `crispr_screen/R/utils.R` is used by gene-list derivation and enrichment filtering; several figure scripts still set local `pval_thresh <- 0.05` (documented in methods; not refactored without byte-identical proof).
- Gene-list consistency can be verified with `Rscript crispr_screen/scripts/00_derive_gene_lists.R --check` (requires MAGeCK summaries); the main figure runner does not itself re-derive lists.
- KEGG annotations can drift if enrichment is recomputed against the live database; committed CSVs are the thesis values of record.
- PDF/DOCX outputs embed per-run timestamps; PNG outputs are the byte-stable check artefacts.
- Metabolomics runner has known non-termination / figure-regeneration quirks (excluded from thesis).
- Exact `renv.lock` for the thesis environment is not reconstructed here.

---

## Citation

Cite this computational repository as supporting material for the Master’s
thesis of **Furkan Emre Bora**, archived at
[https://github.com/femrebora/Thesis](https://github.com/femrebora/Thesis).

Upstream pipeline metadata records the working thesis title
*Metabolic Determinants of TRAIL Resistance in GBM* (Cingöz Lab). Prefer the
title page of the submitted thesis when citing the thesis itself. The thesis
PDF is not distributed in this repository; do not invent a DOI or publication
record.

---

## License / reuse status

All rights reserved. Contact the repository owner for permissions regarding code reuse.

Separately:

- **Laboratory experimental data** are not redistributed through this public repository where restricted.
- **TCGA/GDC data** remain subject to GDC / TCGA access and use terms.
- No open-source license file is granted by this repository at present.
