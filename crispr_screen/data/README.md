# CRISPR screen — input data

The files below are **not distributed** with the repository (they are derived from
restricted sequencing data) and are git-ignored. Place them here to run the
pipeline. The expected schemas are documented so the inputs can be regenerated
from the upstream count/test steps described in `crispr_screen/README.md`.

## `count_table.txt` (tab-separated)

MAGeCK count matrix produced by the `countguides.sh` counting step.

| column | meaning |
|---|---|
| `sgRNA` | guide identifier |
| `Gene`  | target gene symbol (`INTERGENIC`/`CONTROL` for non-targeting) |
| `R0_1..R0_3` | resistant day-0 baseline, biological replicates |
| `Rpost_1..Rpost_4` | resistant tumour-derived samples |
| `S0_1..S0_3` | sensitive day-0 baseline |
| `Spost_1..Spost_3` | sensitive tumour-derived samples |

## `<comparison>.gene_summary.tsv` (tab-separated, MAGeCK `test` output)

One file per comparison: `Rpost_vs_R0`, `Spost_vs_Rpost`, `Spost_vs_S0`.

Standard MAGeCK gene-summary columns:
`id`, `num`, `neg|score`, `neg|p-value`, `neg|fdr`, `neg|rank`, `neg|goodsgrna`,
`neg|lfc`, `pos|score`, `pos|p-value`, `pos|fdr`, `pos|rank`, `pos|goodsgrna`,
`pos|lfc`. The negative tail tests depletion (LFC < 0); the positive tail tests
enrichment (LFC > 0).

## `gene_lists/*.txt`

Plain gene-symbol lists (one per line) used as enrichment inputs. These contain
gene symbols only (no per-patient or raw-count data) and **are** committed.

## Sample-to-group mapping

See `crispr_screen/metadata/sample_metadata_template.csv` for the column layout
used to label the four conditions (S0, Spost, R0, Rpost) and replicates.
