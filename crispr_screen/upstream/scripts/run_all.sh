#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
# run_all.sh — Master Execution Wrapper
# ═══════════════════════════════════════════════════════════════════════════════
# Executes the complete CRISPR-Seq Analysis Pipeline for PRJ25-131.
#
# Usage:
#   bash run_all.sh                    # Full pipeline
#   bash run_all.sh --skip-fastq       # Skip FASTQ preparation
#   bash run_all.sh --skip-count       # Skip MAGeCK count
#   bash run_all.sh --figures-only     # Only regenerate figures
#   bash run_all.sh --new-analysis     # Also run New Analysis (steps 6-8)
#   bash run_all.sh --dry-run          # Print commands without executing
#
# Requirements:
#   conda activate crispr-pipeline
#   (or equivalent Python/R/MAGeCK environment)
# ═══════════════════════════════════════════════════════════════════════════════

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CRISPR_ROOT="${PROJECT_ROOT}/CRISPR-Seq_Analysis"

FASTQ_DIR="${CRISPR_ROOT}/Fastq"
LIB_FILE="${CRISPR_ROOT}/nf-core_crisprseq/library_sabatini_metko.tsv"
COMBINED_DIR="${SCRIPT_DIR}/combined_fastq"
COUNT_DIR="${SCRIPT_DIR}/count"
TEST_DIR="${SCRIPT_DIR}/test"
FIG_DIR="${SCRIPT_DIR}/figures"
REPORT_DIR="${SCRIPT_DIR}/reports"
IMPR_DIR="${SCRIPT_DIR}/improved_results"
SCAF_DIR="${SCRIPT_DIR}/scaffold_results"

PYTHON="python3"
MAGECK="mageck"

# ── Flags ────────────────────────────────────────────────────────────────────
SKIP_FASTQ=false
SKIP_COUNT=false
FIGURES_ONLY=false
NEW_ANALYSIS=false
DRY_RUN=false

for arg in "$@"; do
    case "$arg" in
        --skip-fastq)    SKIP_FASTQ=true ;;
        --skip-count)    SKIP_COUNT=true ;;
        --figures-only)  FIGURES_ONLY=true ;;
        --new-analysis)  NEW_ANALYSIS=true ;;
        --dry-run)       DRY_RUN=true ;;
        --help|-h)
            echo "Usage: bash run_all.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skip-fastq      Skip FASTQ concatenation (step 1)"
            echo "  --skip-count      Skip MAGeCK count (step 2)"
            echo "  --figures-only    Only regenerate figures (step 5)"
            echo "  --new-analysis    Also run New Analysis (steps 6-8)"
            echo "  --dry-run         Print commands without executing"
            echo "  --help, -h        Show this help"
            exit 0
            ;;
    esac
done

# ── Helper ───────────────────────────────────────────────────────────────────
run_step() {
    local step_name="$1"
    local cmd="$2"
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════╗"
    echo "║  ${step_name}"
    echo "╚══════════════════════════════════════════════════════════════════╝"
    echo "  \$ ${cmd}"
    echo ""
    if [ "$DRY_RUN" = false ]; then
        eval "${cmd}" || {
            echo ""
            echo "  [ERROR] Step failed: ${step_name}"
            echo "  Check logs above for details."
            exit 1
        }
    fi
}

# ── Pre-flight Checks ────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════"
echo "  CRISPR-Seq Reproducible Pipeline — PRJ25-131"
echo "  Project: GBM TRAIL Metabolic CRISPR/Cas9 Screen"
echo "  Date:    $(date '+%Y-%m-%d %H:%M:%S')"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Check Python
if ! command -v ${PYTHON} &>/dev/null; then
    echo "[ERROR] Python not found. Activate conda environment first:"
    echo "  conda activate crispr-pipeline"
    exit 1
fi
echo "[OK] Python: $(${PYTHON} --version)"

# Check MAGeCK
if ! command -v ${MAGECK} &>/dev/null; then
    echo "[ERROR] MAGeCK not found. Install: conda install -c bioconda mageck"
    exit 1
fi
echo "[OK] MAGeCK: $(${MAGECK} --version 2>&1 | head -1)"

# Check R
if command -v Rscript &>/dev/null; then
    echo "[OK] R: $(Rscript --version 2>&1)"
else
    echo "[WARN] R not found. R-based figures will be skipped."
fi

echo ""
echo "Flags:"
echo "  SKIP_FASTQ=${SKIP_FASTQ}"
echo "  SKIP_COUNT=${SKIP_COUNT}"
echo "  FIGURES_ONLY=${FIGURES_ONLY}"
echo "  NEW_ANALYSIS=${NEW_ANALYSIS}"
echo "  DRY_RUN=${DRY_RUN}"
echo ""

# ── FIGURES ONLY MODE ────────────────────────────────────────────────────────
if [ "$FIGURES_ONLY" = true ]; then
    run_step "STEP 5 — Generate Figures" \
        "${PYTHON} ${SCRIPT_DIR}/05_generate_figures.py \
            --test-dir ${TEST_DIR} \
            --count-file ${COUNT_DIR}/count_table_fixed.txt \
            --outdir ${FIG_DIR}"
    echo ""
    echo "Figures regenerated: ${FIG_DIR}"
    exit 0
fi

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Prepare FASTQ
# ═══════════════════════════════════════════════════════════════════════════════
if [ "$SKIP_FASTQ" = false ]; then
    run_step "STEP 1 — Concatenate R1+R2 FASTQ" \
        "${PYTHON} ${SCRIPT_DIR}/01_prepare_fastq.py \
            --fastq-dir ${FASTQ_DIR} \
            --outdir ${COMBINED_DIR} \
            --skip-existing"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — MAGeCK Count
# ═══════════════════════════════════════════════════════════════════════════════
if [ "$SKIP_COUNT" = false ]; then
    run_step "STEP 2 — MAGeCK Count" \
        "${PYTHON} ${SCRIPT_DIR}/02_mageck_count.py \
            --fastq-dir ${COMBINED_DIR} \
            --library ${LIB_FILE} \
            --outdir ${COUNT_DIR} \
            --scaffold GTTTTAGAGCTAGAAATAGC \
            --sgrna-len 20"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — MAGeCK Test
# ═══════════════════════════════════════════════════════════════════════════════
run_step "STEP 3 — MAGeCK RRA Test" \
    "${PYTHON} ${SCRIPT_DIR}/03_mageck_test.py \
        --count-file ${COUNT_DIR}/count_table_fixed.txt \
        --outdir ${TEST_DIR}"

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — QC Metrics
# ═══════════════════════════════════════════════════════════════════════════════
run_step "STEP 4 — QC Metrics (Dynamic)" \
    "${PYTHON} ${SCRIPT_DIR}/04_qc_metrics.py \
        --count-file ${COUNT_DIR}/count_table_fixed.txt \
        --outdir ${REPORT_DIR}"

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Generate Figures
# ═══════════════════════════════════════════════════════════════════════════════
run_step "STEP 5 — Generate Publication Figures" \
    "${PYTHON} ${SCRIPT_DIR}/05_generate_figures.py \
        --test-dir ${TEST_DIR} \
        --count-file ${COUNT_DIR}/count_table_fixed.txt \
        --outdir ${FIG_DIR}"

# ═══════════════════════════════════════════════════════════════════════════════
# STEP R — R-based Figures (optional)
# ═══════════════════════════════════════════════════════════════════════════════
if command -v Rscript &>/dev/null; then
    RMD_FILE="${SCRIPT_DIR}/figure_scripts/figures_mageck.Rmd"
    if [ -f "$RMD_FILE" ]; then
        echo ""
        echo "--- R Markdown Figures ---"
        Rscript -e "
            rmarkdown::render('${RMD_FILE}',
                params = list(
                    test_dir = '${TEST_DIR}',
                    count_file = '${COUNT_DIR}/count_table_fixed.txt',
                    output_dir = '${FIG_DIR}'
                ),
                output_dir = '${FIG_DIR}',
                quiet = TRUE
            )
        " || echo "[WARN] R figure generation had issues. Python figures OK."
        echo "[OK] R figures generated"
    fi
fi

# ═══════════════════════════════════════════════════════════════════════════════
# STEPS 6-8 — New Analysis (optional)
# ═══════════════════════════════════════════════════════════════════════════════
if [ "$NEW_ANALYSIS" = true ]; then
    run_step "STEP 6 — Deep QC Inspection" \
        "${PYTHON} ${SCRIPT_DIR}/06_new_analysis_qc.py \
            --mt-dir ${CRISPR_ROOT}/analysis/mageck_test \
            --lt-dir ${CRISPR_ROOT}/analysis/results/last_test \
            --outdir ${FIG_DIR} \
            --report-dir ${REPORT_DIR}"

    run_step "STEP 7 — Improved Analysis" \
        "${PYTHON} ${SCRIPT_DIR}/07_improved_analysis.py \
            --count-file ${CRISPR_ROOT}/analysis/results/last_count/last_count_fixed.txt \
            --outdir ${IMPR_DIR} \
            --fig-dir ${FIG_DIR}"

    run_step "STEP 8 — Scaffold Comparison" \
        "${PYTHON} ${SCRIPT_DIR}/08_scaffold_comparison.py \
            --fastq-dir ${FASTQ_DIR} \
            --outdir ${SCAF_DIR} \
            --fig-dir ${FIG_DIR}"
fi

# ═══════════════════════════════════════════════════════════════════════════════
# DONE
# ═══════════════════════════════════════════════════════════════════════════════
echo ""
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║  PIPELINE COMPLETE                                              ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
echo "Output directories:"
echo "  Combined FASTQ:   ${COMBINED_DIR}"
echo "  MAGeCK count:     ${COUNT_DIR}"
echo "  MAGeCK test:      ${TEST_DIR}"
echo "  QC reports:       ${REPORT_DIR}"
echo "  Figures:          ${FIG_DIR}"
if [ "$NEW_ANALYSIS" = true ]; then
    echo "  Improved results: ${IMPR_DIR}"
    echo "  Scaffold results: ${SCAF_DIR}"
fi
echo ""
echo "Figure count: $(ls ${FIG_DIR}/*.pdf 2>/dev/null | wc -l) PDFs, $(ls ${FIG_DIR}/*.png 2>/dev/null | wc -l) PNGs"
echo ""
echo "Next steps:"
echo "  1. Review QC report:    less ${REPORT_DIR}/qc_report.txt"
echo "  2. Check QC metrics:    less ${REPORT_DIR}/qc_metrics.csv"
echo "  3. Review figures:      ls ${FIG_DIR}/"
echo "  4. Check significant:   less ${FIG_DIR}/significant_genes.csv"
echo ""
echo "For biological interpretation, use:"
echo "  - string_database skill for PPI networks"
echo "  - reactome_database skill for pathway enrichment"
echo "  - uniprot_database skill for protein function"
echo ""
echo "Date completed: $(date '+%Y-%m-%d %H:%M:%S')"
