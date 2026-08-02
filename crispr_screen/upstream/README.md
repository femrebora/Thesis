# crispr_screen/upstream/ — FASTQ → count matrix → MAGeCK gene summaries

Everything that runs **before** the R figure scripts in `crispr_screen/scripts/`.
The downstream figure pipeline consumes only the three MAGeCK gene-summary files
and the count matrix that this stage produces.

```
FASTQ (paired-end, 13 samples)
   │
   ├─ nf-core/CRISPR_nfcore.sh ......... nf-core/crisprseq v2.3.0 (QC + first-pass count/test)
   │
   └─ scripts/run_all.sh ............... the pipeline of record
        01_prepare_fastq.py ........... concatenate R1+R2 per sample
        02_mageck_count.py ............ MAGeCK count vs Sabatini metabolic library
        03_mageck_test.py ............. MAGeCK RRA test, 3 comparisons
        04_qc_metrics.py .............. Gini / zero-fraction / detection rate
        05_generate_figures.py ........ first-generation Python figures (superseded)
        06_new_analysis_qc.py ......... nf-core (R1-only) vs R1+R2 comparison
        07_improved_analysis.py ....... parameter-sensitivity re-runs
        08_scaffold_comparison.py ..... 20 bp vs 19 bp scaffold effect
                    │
                    ▼
        count_table.txt + <comparison>.gene_summary.tsv
                    │
                    ▼
        crispr_screen/scripts/*.R  (thesis figures)
```

`REPRODUCIBLE_PIPELINE.md` is the full protocol: sample groups, comparison
definitions, environment setup, QC thresholds and known caveats.

## Two counting routes, one of record

`nf-core/` ran crisprseq v2.3.0 and counted **R1 only**. That undercounts,
because ~34.6% of sgRNAs are sequenced in reverse-complement orientation on R1.
The pipeline of record therefore concatenates R1+R2 (`01_prepare_fastq.py`) and
re-counts with `02_mageck_count.py`. `06_new_analysis_qc.py` quantifies the
difference between the two routes.

**The gene summaries consumed by the thesis figures come from the R1+R2 route.**

## Relationship to the figures in the thesis

`05_generate_figures.py` produced the first generation of figures in Python.
Those were later rebuilt in R (`crispr_screen/scripts/01–12`) for the Turkish
thesis with consistent typography and a shared palette. **The R scripts are the
source of the submitted figures**; the Python version is kept here as the
provenance record, not as a parallel figure pipeline.

## Inputs not distributed here

FASTQ files, the sgRNA library TSV, `samplesheet.csv` and `contrasts.csv` are
derived from restricted sequencing data and are git-ignored. See
`crispr_screen/data/README.md` for the expected schemas.

## Reproducing this stage

```bash
conda activate crispr-pipeline      # needs mageck >= 0.5.9.5, python >= 3.10
cd crispr_screen/upstream/scripts
bash run_all.sh --dry-run           # print the commands first
bash run_all.sh
```

Requires the restricted inputs to be present. If you only need to regenerate
figures, skip this stage entirely — the MAGeCK outputs it produces are the
documented entry point for `crispr_screen/scripts/run_crispr_analysis.R`.
