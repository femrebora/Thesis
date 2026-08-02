# CRISPR-Seq Reproducible Analysis Pipeline

**Project:** PRJ25-131 — GBM TRAIL Metabolic CRISPR/Cas9 Screen  
**Library:** Sabatini Human Metabolic Gene Knockout Library (Addgene #110066)  
**Lab:** Cingöz Lab  
**Version:** v1.0.0 — 2026-06-09  
**Protocol Reference:** THESIS_MASTER_PROTOCOL_v5.13_FINAL.md

---

## Table of Contents

1. [Overview](#overview)
2. [Skills, Plugins & MCP Routing](#skills-plugins--mcp-routing)
3. [Environment Setup](#environment-setup)
4. [Input Data](#input-data)
5. [Pipeline Architecture](#pipeline-architecture)
6. [Step-by-Step Execution](#step-by-step-execution)
7. [Figure Inventory](#figure-inventory)
8. [QC Thresholds & Reporting](#qc-thresholds--reporting)
9. [Known Issues & Caveats](#known-issues--caveats)
10. [Final Submission Checklist](#final-submission-checklist)

---

## Overview

This pipeline performs end-to-end analysis of a pooled CRISPR/Cas9 metabolic gene knockout screen. It processes paired-end FASTQ files from 13 samples across 4 experimental conditions (S0, Spost, R0, Rpost), runs MAGeCK for sgRNA counting and statistical testing, computes QC metrics from count data (never hardcoded), and generates publication-quality figures.

### Pipeline Summary

```
FASTQ (PE) → Concat R1+R2 → MAGeCK count → MAGeCK test (3 comparisons)
                                         → QC metrics (dynamic)
                                         → Figures (9 main + QC + New Analysis)
```

### Sample Groups

| Condition | Samples | Description |
|-----------|---------|-------------|
| S0 | S01, S02, S03 | Sensitive cells — Day 0 baseline |
| Spost | S1, S2, S3 | Sensitive cells — Post-TRAIL treatment |
| R0 | R01, R02, R03 | Resistant cells — Day 0 baseline |
| Rpost | R1, R2, R3, R4 | Resistant cells — Post-TRAIL treatment |

### Comparisons

| Comparison | Label | Treatment | Control |
|-----------|-------|-----------|---------|
| Spost_vs_S0 | Sensitive: Post vs Day-0 | S1,S2,S3 | S01,S02,S03 |
| Rpost_vs_R0 | Resistant: Post vs Day-0 | R1,R2,R3,R4 | R01,R02,R03 |
| Spost_vs_Rpost | Sensitive vs Resistant (Post) | S1,S2,S3 | R1,R2,R3,R4 |

---

## Skills, Plugins & MCP Routing

This pipeline integrates with the Claude Code skills/plugins ecosystem for thesis work. When invoked, the following are available:

### Active MCP Servers

| MCP Server | Purpose | Pipeline Phase |
|-----------|---------|---------------|
| `context-mode` | Context-efficient file ops, sandboxed code execution | All phases |
| `elicit` | Literature search for gene function annotation | Phase 4 (interpretation) |
| `claude-mem` | Session memory, observation search | All phases |

### Available Skills

| Skill | Purpose | When to Invoke |
|-------|---------|---------------|
| `bioinformatics` | Variant scoring, pipeline standards, QC thresholds | QC assessment, pipeline design |
| `avoid-ai-writing` | AI-writing detection avoidance | Report generation |
| `opendraft-drafting` | Thesis drafting prompts (24 templates) | Writing results interpretation |
| `academic-writing-agents:academic` | Multi-agent academic writing pipeline | Full paper drafting |
| `literature_search_europepmc` | Europe PMC paper search | Gene function lookups |
| `pubmed_database` | PubMed querying | Literature context |
| `string_database` | Protein-protein interaction networks | Pathway analysis of hits |
| `reactome_database` | Pathway database | Hit gene pathway enrichment |
| `uniprot_database` | Protein function annotation | Hit gene characterization |

### Plugin Agents (academic-writing-agents v3.7.0)

| Agent | Role | Pipeline Phase |
|-------|------|---------------|
| `synthesis_agent` | Cross-source integration, contradiction resolution | Phase 5 (report) |
| `research_architect_agent` | Methodology blueprint | Phase 0 (planning) |
| `report_compiler_agent` | APA 7.0 report drafting | Phase 5 (report) |
| `bibliography_agent` | Reference management | Phase 5 (report) |
| `prose-polisher` | Academic prose improvement | Phase 5 (report) |
| `technical-reviewer` | Technical accuracy verification | Phase 6 (review) |
| `consistency-checker` | Terminology/cross-reference check | Phase 6 (review) |
| `logic-reviewer` | Logical flow, argument structure | Phase 6 (review) |

### Tool Routing Matrix

```
PHASE 1 — FASTQ Prep      → Bash, Python (no agents needed)
PHASE 2 — MAGeCK count     → Bash, Python (no agents needed)
PHASE 3 — MAGeCK test      → Bash, Python, R (no agents needed)
PHASE 4 — QC & Figures     → Python, R, bioinformatics skill
PHASE 5 — Interpretation   → research-analyst, synthesis_agent, string_database, reactome_database
PHASE 6 — Report Writing   → report_compiler_agent, prose-polisher, technical-reviewer
PHASE 7 — Review & Submit  → consistency-checker, logic-reviewer, avoid-ai-writing
```

---

## Environment Setup

### Option A: Conda (Recommended)

```bash
# Create environment
conda create -n crispr-pipeline -c bioconda -c conda-forge \
    python=3.11 \
    mageck=0.5.9.5 \
    numpy pandas matplotlib seaborn scipy \
    adjusttext \
    r-base=4.3 \
    r-rmarkdown r-knitr r-ggplot2 r-dplyr r-tidyr \
    r-plotly r-pheatmap r-ggrepel \
    fastqc multiqc \
    -y

conda activate crispr-pipeline

# Python packages not in conda
pip install adjustText

# Verify installation
mageck --version
python -c "import matplotlib; print('OK')"
Rscript -e "library(ggplot2); print('OK')"
```

### Option B: Docker

```bash
# Build image
docker build -t crispr-pipeline -f Dockerfile .

# Run pipeline
docker run --rm -v $(pwd):/data crispr-pipeline bash run_all.sh
```

### Option C: Singularity/Apptainer (HPC)

```bash
apptainer build crispr-pipeline.sif docker://emrebora/crispr-pipeline:latest
apptainer run --bind /data crispr-pipeline.sif bash run_all.sh
```

### Required External Tools

| Tool | Version | Purpose |
|------|---------|---------|
| MAGeCK | ≥0.5.9.5 | sgRNA counting and statistical testing |
| Python | ≥3.10 | Pipeline orchestration, figure generation |
| R | ≥4.3 | MAGeCK RRA visualization, QC reports |
| FastQC | ≥0.12 | Raw read quality control |
| MultiQC | ≥1.21 | QC report aggregation |
| Nextflow | ≥22.10 | nf-core/crisprseq pipeline execution |

---

## Input Data

### Required Files

```
CRISPR-Seq_Analysis/
├── Fastq/
│   ├── PRJ25-131-R01_L01_59_1.fq.gz    # R01 R1
│   ├── PRJ25-131-R01_L01_59_2.fq.gz    # R01 R2
│   ├── ... (26 FASTQ files total; 13 samples × 2 reads)
│   └── PRJ25-131-S3_L01_64_2.fq.gz     # S3 R2
├── nf-core_crisprseq/
│   └── library_sabatini_metko.tsv       # sgRNA library (sgRNA, sequence, gene)
```

### Library Format (`library_sabatini_metko.tsv`)

TSV with columns: `sgRNA_ID`, `sgRNA_sequence`, `Gene`

```
sgRNA_ID    sgRNA_sequence          Gene
A1BG_sg1    ACGT...                 A1BG
A1BG_sg2    TGCA...                 A1BG
...
```

### FASTQ Naming Convention

```
PRJ25-131-{sample_id}_L01_{numeric_id}_{read}.fq.gz
```

### Sample Sheet (for nf-core pipeline)

```csv
sample,fastq_1,fastq_2,condition,trim_5
PRJ25-131-R01_L01_59,Fastq/PRJ25-131-R01_L01_59_1.fq.gz,Fastq/PRJ25-131-R01_L01_59_2.fq.gz,R0,GTGGAAAGGACGAAACACCG
...
```

The `trim_5` column contains the U6 promoter end sequence (GTGGAAAGGACGAAACACCG) for lentiCRISPR v1, used to locate sgRNA sequences in reads.

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                     CRISPR-Seq Analysis Pipeline                     │
│                     PRJ25-131 GBM TRAIL Screen                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  ┌──────────┐   ┌──────────────┐   ┌──────────────┐                │
│  │  FASTQ   │──▶│ 01_prepare   │──▶│ 02_mageck    │                │
│  │  (PE)    │   │ _fastq.py    │   │ _count.py    │                │
│  │  26 files│   │ Concat R1+R2 │   │ sgRNA counts │                │
│  └──────────┘   └──────────────┘   └──────┬───────┘                │
│                                            │                         │
│                    ┌───────────────────────┤                         │
│                    ▼                       ▼                         │
│            ┌──────────────┐      ┌──────────────────┐              │
│            │ 03_mageck    │      │ 04_qc_metrics.py │              │
│            │ _test.py     │      │ Dynamic QC from   │              │
│            │ 3 comparisons│      │ count data        │              │
│            └──────┬───────┘      └────────┬─────────┘              │
│                   │                       │                          │
│                   └───────────┬───────────┘                          │
│                               ▼                                      │
│                     ┌──────────────────┐                            │
│                     │ 05_generate      │                            │
│                     │ _figures.py      │                            │
│                     │ 9 main figures   │                            │
│                     └──────────────────┘                            │
│                               │                                      │
│            ┌──────────────────┼──────────────────┐                  │
│            ▼                  ▼                  ▼                   │
│   ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐        │
│   │ 06_new_      │  │ 07_improved  │  │ 08_scaffold      │        │
│   │ analysis_qc  │  │ _analysis.py │  │ _comparison.py   │        │
│   │ Deep QC      │  │ Better params│  │ 19bp vs 20bp     │        │
│   └──────────────┘  └──────────────┘  └──────────────────┘        │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Step-by-Step Execution

### Quick Start (Full Pipeline)

```bash
conda activate crispr-pipeline
cd crispr_screen/upstream/scripts
bash run_all.sh
```

### Individual Steps

#### Step 0 — Verify Environment

```bash
python -c "
import numpy, pandas, matplotlib, seaborn, scipy
print('Python environment OK')
"
mageck --version
Rscript -e 'cat(paste0(\"R \", R.version\$major, \".\", R.version\$minor, \"\\n\"))'
```

#### Step 1 — Prepare FASTQ (Concatenate R1+R2)

**Script:** `01_prepare_fastq.py`  
**Input:** `CRISPR-Seq_Analysis/Fastq/*.fq.gz`  
**Output:** `results/combined_fastq/{label}_combined.fq.gz`

```bash
python 01_prepare_fastq.py \
    --fastq-dir ../CRISPR-Seq_Analysis/Fastq \
    --outdir combined_fastq
```

**Rationale:** ~34.6% of R1 reads carry sgRNA in reverse-complement orientation. R2 reads have forward orientation. Concatenating both recovers sgRNAs in both orientations, improving detection rate.

#### Step 2 — MAGeCK Count

**Script:** `02_mageck_count.py`  
**Input:** Combined FASTQ files + library TSV  
**Output:** `results/count/count_table.txt`, `results/count/count_table_fixed.txt`

```bash
python 02_mageck_count.py \
    --fastq-dir combined_fastq \
    --library ../CRISPR-Seq_Analysis/nf-core_crisprseq/library_sabatini_metko.tsv \
    --outdir count \
    --scaffold GTTTTAGAGCTAGAAATAGC \
    --sgrna-len 20
```

**Parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--scaffold` | GTTTTAGAGCTAGAAATAGC | 20bp lentiCRISPRv1 scaffold (includes leading G) |
| `--sgrna-len` | 20 | sgRNA length (standard for GeCKO/Sabatini libraries) |
| `--norm-method` | median | Normalization method (median preferred for skewed data) |

#### Step 3 — MAGeCK Test (RRA)

**Script:** `03_mageck_test.py`  
**Input:** Fixed count table  
**Output:** `results/test/{comparison}.gene_summary.txt`, `results/test/{comparison}.sgrna_summary.txt`

```bash
python 03_mageck_test.py \
    --count-file count/count_table_fixed.txt \
    --outdir test
```

**Comparisons performed:**

1. `Spost_vs_S0`: Sensitive post-treatment vs Day-0 (S1,S2,S3 vs S01,S02,S03)
2. `Rpost_vs_R0`: Resistant post-treatment vs Day-0 (R1,R2,R3,R4 vs R01,R02,R03)
3. `Spost_vs_Rpost`: Sensitive vs Resistant (S1,S2,S3 vs R1,R2,R3,R4)

**MAGeCK test parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--gene-lfc-method` | median | Per-gene LFC from median sgRNA LFC |
| `--remove-zero` | both | Remove sgRNAs with zero counts in both groups |
| `--remove-zero-threshold` | 0 | Threshold for zero-count removal |

**Important:** This pipeline uses `--remove-zero both` with threshold 0. For improved sensitivity with more testable genes, use `07_improved_analysis.py` which uses `--remove-zero-threshold 10`.

#### Step 4 — QC Metrics (Dynamically Computed)

**Script:** `04_qc_metrics.py`  
**Input:** Count table + gene summaries  
**Output:** `results/reports/qc_metrics.csv`

```bash
python 04_qc_metrics.py \
    --count-file count/count_table_fixed.txt \
    --test-dir test \
    --outdir reports
```

**QC metrics computed dynamically (NEVER hardcoded):**

| Metric | Formula | Ideal | Poor |
|--------|---------|-------|------|
| Mapping rate (%) | Mapped / Total × 100 | ≥60% | <10% |
| Gini index | Σ(2i-n-1)xᵢ / (n·Σxᵢ) | <0.2 | >0.9 |
| Zero-count sgRNAs (%) | N(zero) / N(total) × 100 | <20% | >80% |
| Detected sgRNAs | N(nonzero) | >27,000 (90%) | <3,000 (10%) |
| Reads per sample | Σ counts | >9M | <3M |

#### Step 5 — Generate Publication Figures

**Script:** `05_generate_figures.py`  
**Input:** Gene summaries + count table  
**Output:** `results/figures/Fig{1-9}_*.{pdf,png}`

```bash
python 05_generate_figures.py \
    --test-dir test \
    --count-file count/count_table_fixed.txt \
    --outdir figures \
    --pval-threshold 0.05 \
    --top-label 10
```

**Figure list (see [Figure Inventory](#figure-inventory))**

#### Step 6 — Deep QC Inspection (New Analysis)

**Script:** `06_new_analysis_qc.py`  
**Input:** Original + combined run results  
**Output:** `results/figures/QC{1-5}_*.{pdf,png}`, `results/reports/qc_comparison_report.txt`

```bash
python 06_new_analysis_qc.py \
    --mt-dir ../CRISPR-Seq_Analysis/analysis/mageck_test \
    --lt-dir ../CRISPR-Seq_Analysis/analysis/results/last_test \
    --outdir figures \
    --report-dir reports
```

#### Step 7 — Improved Analysis

**Script:** `07_improved_analysis.py`  
**Input:** Count tables  
**Output:** `results/improved_results/`, `results/figures/IMP{1,2}_*.{pdf,png}`

```bash
python 07_improved_analysis.py \
    --count-file ../CRISPR-Seq_Analysis/analysis/results/last_count/last_count_fixed.txt \
    --outdir improved_results \
    --fig-dir figures
```

**Improvements:**
1. `--remove-zero-threshold 10` (more testable genes)
2. `--norm-method none` (avoid double normalization)
3. Two-sided p-value correction (×2 Bonferroni)
4. Analysis with and without Spost_1 outlier

#### Step 8 — Scaffold Comparison

**Script:** `08_scaffold_comparison.py`  
**Input:** FASTQ files  
**Output:** `results/figures/SCAF1_*.{pdf,png}`, `results/scaffold_results/`

```bash
python 08_scaffold_comparison.py \
    --fastq-dir ../CRISPR-Seq_Analysis/Fastq \
    --outdir scaffold_results \
    --fig-dir figures
```

---

## Figure Inventory

### Main Figures (05_generate_figures.py)

| # | Name | Description | Type |
|---|------|-------------|------|
| 1 | `Fig1_volcano` | Volcano plots (LFC vs −log₁₀p) — 3 panels | Scatter |
| 2 | `Fig2_rank_lfc` | Gene rank by fold change — 3 panels | Scatter |
| 3 | `Fig3_heatmap` | LFC heatmap across all comparisons | Heatmap |
| 4 | `Fig4_barplot` | Top 15 significant genes per comparison | Bar |
| 5 | `Fig5_qc_metrics` | QC panels: depth, zero-count, Gini, detected | Multi-panel bar |
| 6 | `Fig6_sgrna_profiles` | Per-sgRNA RPM profiles for key hit genes | Line |
| 7 | `Fig7_replicate_corr` | Replicate Spearman correlations | Scatter matrix |
| 8 | `Fig8_cross_comparison` | Bubble plot across comparisons | Bubble |
| 9 | `Fig9_summary_table` | Table of all significant genes with annotations | Table |

### New Analysis QC Figures (06_new_analysis_qc.py)

| # | Name | Description |
|---|------|-------------|
| QC1 | `QC1_mapping_rate_comparison` | nf-core vs R1+R2 mapping rates |
| QC2 | `QC2_gini_index` | Gini index per sample (combined run) |
| QC3 | `QC3_testable_genes_comparison` | Testable gene count comparison |
| QC4 | `QC4_lfc_comparison_{comp}` | LFC & p-value scatter between runs (×3) |
| QC5 | `QC5_nfcore_count_dist` / `QC5_combined_count_dist` | sgRNA count distributions |

### Improved Analysis Figures (07_improved_analysis.py)

| # | Name | Description |
|---|------|-------------|
| IMP1 | `IMP1_sgrna_profiles_hits` | Per-sgRNA RPM for 8 key hit genes |
| IMP2 | `IMP2_volcano_improved` | Volcano with corrected p-values |

### Scaffold Comparison Figures (08_scaffold_comparison.py)

| # | Name | Description |
|---|------|-------------|
| SCAF1 | `SCAF1_scaffold_comparison` | 20bp vs 19bp scaffold detection rates |

---

## QC Thresholds & Reporting

### QC Grade Classification

| Grade | Criteria | Action |
|-------|----------|--------|
| **A — Pass** | All metrics in ideal range | Proceed to biological interpretation |
| **B — Acceptable** | 1-2 metrics in warning range | Note as limitations; use p < 0.05 |
| **C — Caution** | >2 metrics in warning range | Flag as exploratory; require orthogonal validation |
| **D — Fail** | >2 metrics critical | Do NOT draw conclusions; experiment repeat required |

### Current Screen Assessment (PRJ25-131)

| Metric | Value | Ideal | Grade |
|--------|-------|-------|-------|
| Mapping rate | 0.4–6.5% | >60% | 🔴 **D** |
| Gini index | 0.73–0.99 | <0.2 | 🔴 **D** |
| Zero-sgRNA % | >99% (all samples) | <20% | 🔴 **D** |
| Detected sgRNAs | 1,006 / 30,197 (3.3%) | >90% | 🔴 **D** |
| Read depth | 9.5K–28.3M | >9M | 🟡 B |
| Duplication rate | 89–97% | <30% | 🔴 **D** |

**Overall Grade: D — Data quality is insufficient for definitive biological conclusions. All results are exploratory.**

### QC Reporting Template

```markdown
## QC Report — PRJ25-131 — {date}

### Sample-level metrics
| Sample | Condition | Total Reads (M) | Mapped (M) | Mapping % | Gini | Zero % | Detected sgRNAs |
|--------|-----------|-----------------|------------|-----------|------|--------|-----------------|
| ... | ... | ... | ... | ... | ... | ... | ... |

### Comparison-level metrics
| Comparison | Testable Genes | p<0.05 | Min FDR | Max |LFC| |
|------------|---------------|--------|---------|----------|
| ... | ... | ... | ... | ... |

### Flags
- [ ] Mapping rate <10% on N samples
- [ ] Gini >0.9 on N samples
- [ ] S1 excluded (only 5 sgRNAs detected)
```

---

## Known Issues & Caveats

### 1. Severe Library Bottleneck (CRITICAL)

The sequencing library suffered from extreme bottlenecking during the wet-lab phase:
- Only 3.3% of sgRNAs were detected (1,006 out of 30,197)
- Mapping rates 0.4–6.5% (ideal: >60%)
- Gini indices 0.73–0.99 (ideal: <0.2)
- Duplication rates 89–97%

**Impact:** This is a wet-lab issue that cannot be fixed bioinformatically. The experiment should be repeated with:
- ≥15 million cells at transduction (500× coverage)
- MOI 0.3–0.5
- lentiCRISPR v2 (Addgene #52961)
- Plasmid library NGS QC before viral production

### 2. S1 Sample Unusable

Spost_1 (S1) has only 5 non-zero sgRNAs (Gini = 1.000). This sample is excluded from improved analyses.

### 3. FDR Cannot Be Applied

Due to the severe bottleneck, FDR values are uniformly high (>0.25 for all genes). The pipeline uses raw p < 0.05 as the primary criterion, but these are exploratory — not confirmatory.

### 4. Single-sgRNA Limitations

Many significant hits are driven by a single sgRNA (where other sgRNAs targeting the same gene had zero counts). Single-sgRNA effects cannot be confirmed as gene-specific. All such genes are flagged with ⚠ in figures.

### 5. Hardcoded QC (HISTORICAL — FIXED)

Previous scripts (`last.py`, `03_downstream_analysis.py`) contained hardcoded QC values. The current pipeline (`pipeline.py` and all scripts in this `results/` directory) computes ALL QC metrics dynamically from count data.

### 6. LFC Direction Bug (HISTORICAL — FIXED)

Previous `downstream_analysis.py` compared `neg|lfc < pos|lfc` to determine direction — these columns always contain identical values, causing all genes to be classified as "enriched." Fixed by using `LFC < 0 → Depleted, LFC > 0 → Enriched`.

### 7. Duplicate Column Names in nf-core Output

nf-core MAGeCK count labels all sample columns identically. `02_mageck_count.py` automatically renames duplicates.

---

## Final Submission Checklist

Before submitting results, verify:

### Reproducibility
- [ ] All scripts run without errors from clean environment
- [ ] Environment YAML specifies exact versions (no `latest`)
- [ ] Random seeds are set where applicable
- [ ] Input data paths are documented
- [ ] All output file paths are documented

### QC
- [ ] QC metrics computed dynamically (verify: no hardcoded values in source)
- [ ] QC report generated and reviewed
- [ ] S1 sample flagged/excluded
- [ ] All single-sgRNA hits flagged
- [ ] Mapping rate, Gini, zero-count % reported per sample

### Figures
- [ ] All figures saved as both PDF and PNG at 300 dpi
- [ ] Figure captions include sample sizes and statistical thresholds
- [ ] Color schemes are colorblind-friendly (DejaVu Sans, Blue/Red/Grey palette)
- [ ] FDR note included where applicable

### Statistics
- [ ] p-value threshold clearly stated
- [ ] Multiple testing correction described (or rationale for omission)
- [ ] Effect sizes (LFC) reported alongside p-values
- [ ] Number of testable genes stated per comparison

### Interpretation
- [ ] All significant genes annotated with known functions
- [ ] Single-sgRNA vs multi-sgRNA distinction made
- [ ] Exploratory nature of results clearly stated
- [ ] Orthogonal validation recommended for top hits

### Code Quality
- [ ] No hardcoded values (all parameters at top of scripts or CLI args)
- [ ] Functions have docstrings
- [ ] Error handling for missing files
- [ ] Progress messages printed during execution

### Documentation
- [ ] README explains full workflow
- [ ] Input/output specifications documented
- [ ] Known issues documented
- [ ] Environment setup instructions complete

---

## References

1. Li W, et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biology*, 15, 554.
2. Li W, et al. (2015). Quality control, modeling, and visualization of CRISPR screens with MAGeCK-VISPR. *Genome Biology*, 16, 281.
3. Wang B, et al. (2019). Integrative analysis of pooled CRISPR genetic screens using MAGeCKFlute. *Nature Protocols*, 14, 756–780.
4. Sanson KR, et al. (2018). Optimized libraries for CRISPR-Cas9 genetic screens with multiple modalities. *Nature Communications*, 9, 5416.
5. nf-core/crisprseq: https://nf-co.re/crisprseq
6. Sabatini Metabolic Library: Addgene #110066 — https://www.addgene.org/110066/

---

## Citation

If you use this pipeline, please cite:

```
Emre Bora. (2026). CRISPR-Seq Reproducible Analysis Pipeline v1.0.0.
PRJ25-131 — GBM TRAIL Metabolic CRISPR/Cas9 Screen.
Cingöz Lab. https://github.com/emrebora/Thesis
```

---

*Pipeline v1.0.0 — Generated 2026-06-09*  
*Protocol compatible with THESIS_MASTER_PROTOCOL_v5.13_FINAL.md*  
*All QC metrics dynamically computed — zero hardcoded values*
