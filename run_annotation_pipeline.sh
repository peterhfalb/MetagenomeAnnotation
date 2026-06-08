#!/bin/bash
# =============================================================================
# run_annotation_pipeline.sh
# Master submission script: set your paths once, submit the full pipeline.
#
# Submits all steps in dependency order. Steps 3-7 run in parallel after
# MetaEuk (step 2) completes. Database setup (step 0) is skipped automatically
# if databases are already built.
#
# Usage:
#   bash run_annotation_pipeline.sh              # normal run
#   bash run_annotation_pipeline.sh --test 2     # test on first 2 samples
#   bash run_annotation_pipeline.sh --force-db   # force re-run database setup
#
# If a step fails mid-run, see the RECOVERY section at the bottom of this file.
#
# After all jobs complete, run step 8 manually (Rscript 08_integrate.R).
# The exact command is printed at submission time.
# =============================================================================

set -euo pipefail

# =============================================================================
# !! EDIT THESE TWO PATHS BEFORE RUNNING !!
# =============================================================================

# Output directory from the MetagenomeAssembly pipeline
ASM_DIR="/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Assembled"

# Where annotation outputs will be written (created automatically)
ANN_DIR="/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation"

# =============================================================================
# Configuration — leave as-is unless you have a reason to change
# =============================================================================

DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_N=""
FORCE_DB=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test)     TEST_N="${2:-2}"; shift 2 ;;
        --force-db) FORCE_DB=true;   shift   ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

TEST_FLAG=""
[[ -n "${TEST_N}" ]] && TEST_FLAG="--test ${TEST_N}"

# ── Validate paths ────────────────────────────────────────────────────────────
if [[ "${ASM_DIR}" == "/path/to/MetagenomeAssembly/output" || \
      "${ANN_DIR}" == "/path/to/annotation/output" ]]; then
    echo "ERROR: Edit ASM_DIR and ANN_DIR at the top of this script before running." >&2
    exit 1
fi
if [[ ! -d "${ASM_DIR}" ]]; then
    echo "ERROR: Assembly directory not found: ${ASM_DIR}" >&2
    exit 1
fi

echo "============================================================"
echo "Annotation pipeline submission"
echo "Assembly dir    : ${ASM_DIR}"
echo "Annotation dir  : ${ANN_DIR}"
[[ -n "${TEST_N}" ]] && echo "TEST MODE       : first ${TEST_N} samples"
echo "============================================================"
echo ""

# =============================================================================
# STEP 0: Database setup — auto-skipped if already complete
# =============================================================================

DB_FLAG="${DB_DIR}/databases_complete.flag"
STEP0_ID=""

check_databases() {
    local all_present=true
    echo "  Database status:"
    local -A checks=(
        ["DIAMOND NR"]="${DB_DIR}/diamond/nr.dmnd"
        ["NCBI taxonomy"]="${DB_DIR}/taxonomy/nodes.dmp"
        ["KOfam profiles"]="${DB_DIR}/kofam/profiles"
        ["dbCAN3"]="${DB_DIR}/dbcan/dbCAN.hmm.h3i"
        ["MetaEuk target DB"]="${DB_DIR}/metaeuk/fungi_refseq_db"
    )
    # Preserve insertion order via explicit list
    for label in "DIAMOND NR" "NCBI taxonomy" "KOfam profiles" "dbCAN3" "MetaEuk target DB"; do
        local path="${checks[$label]}"
        if [[ -e "${path}" ]]; then
            echo "    [OK]     ${label}"
        else
            echo "    [MISSING] ${label}"
            all_present=false
        fi
    done
    # PHI-base checked separately (manual download)
    if [[ -f "${DB_DIR}/phibase/phi-base.dmnd" ]]; then
        echo "    [OK]     PHI-base"
    else
        echo "    [MISSING] PHI-base  (manual download required — see 00_setup_databases.sh Section 5)"
    fi
    ${all_present}
}

if [[ "${FORCE_DB}" == true ]]; then
    echo "Step 0: --force-db set — submitting database setup regardless of current state."
    check_databases || true
    echo ""
    STEP0_ID=$(bash "${SCRIPT_DIR}/submit_setup_databases.sh" | grep -oP '(?<=job: )\d+')
    echo "  Step 0 job ID: ${STEP0_ID}"
    echo ""
elif [[ -f "${DB_FLAG}" ]]; then
    echo "Step 0: databases already complete (found ${DB_FLAG}) — skipping."
    echo ""
else
    echo "Step 0: databases_complete.flag not found — checking individual databases..."
    if check_databases; then
        echo ""
        echo "  All databases present but flag missing. Touching flag and skipping setup."
        touch "${DB_FLAG}"
    else
        echo ""
        echo "  One or more databases missing — submitting database setup job."
        echo "  NOTE: PHI-base requires a manual download first (see 00_setup_databases.sh Section 5)."
        echo "  If PHI-base is the only missing database, touch ${DB_FLAG} and re-run with --force-db=false."
        echo ""
        STEP0_ID=$(bash "${SCRIPT_DIR}/submit_setup_databases.sh" | grep -oP '(?<=job: )\d+')
        echo "  Step 0 job ID: ${STEP0_ID}"
    fi
    echo ""
fi

# =============================================================================
# STEP 1: Tiara contig classification
# =============================================================================

echo "Submitting Step 1: Tiara contig classification..."
STEP1_ARGS=(
    --assembly-dir   "${ASM_DIR}"
    --annotation-dir "${ANN_DIR}"
)
[[ -n "${TEST_FLAG}" ]] && STEP1_ARGS+=( ${TEST_FLAG} )
[[ -n "${STEP0_ID}"  ]] && STEP1_ARGS+=( --after "${STEP0_ID}" )

STEP1_ID=$(bash "${SCRIPT_DIR}/submit_classify_contigs.sh" "${STEP1_ARGS[@]}" \
    | grep -oP '(?<=job: )\d+')
echo "  Job ID: ${STEP1_ID}"
echo ""

# =============================================================================
# STEP 2: MetaEuk gene prediction
# =============================================================================

echo "Submitting Step 2: MetaEuk gene prediction..."
STEP2_ARGS=(
    --assembly-dir   "${ASM_DIR}"
    --annotation-dir "${ANN_DIR}"
    --after          "${STEP1_ID}"
)
[[ -n "${TEST_FLAG}" ]] && STEP2_ARGS+=( ${TEST_FLAG} )

STEP2_ID=$(bash "${SCRIPT_DIR}/submit_metaeuk.sh" "${STEP2_ARGS[@]}" \
    | grep -oP '(?<=job: )\d+')
echo "  Job ID: ${STEP2_ID}"
echo ""

# =============================================================================
# STEPS 3-7: Functional annotation + quantification (all parallel after step 2)
# =============================================================================

echo "Submitting Steps 3-7 (parallel, all depend on Step 2: ${STEP2_ID})..."

PARALLEL_ARGS=( --annotation-dir "${ANN_DIR}" --after "${STEP2_ID}" )
[[ -n "${TEST_FLAG}" ]] && PARALLEL_ARGS+=( ${TEST_FLAG} )

STEP3_ID=$(bash "${SCRIPT_DIR}/submit_kofam.sh"      "${PARALLEL_ARGS[@]}" | grep -oP '(?<=job: )\d+')
STEP4_ID=$(bash "${SCRIPT_DIR}/submit_dbcan.sh"      "${PARALLEL_ARGS[@]}" | grep -oP '(?<=job: )\d+')
STEP5_ID=$(bash "${SCRIPT_DIR}/submit_phibase.sh"    "${PARALLEL_ARGS[@]}" | grep -oP '(?<=job: )\d+')
STEP6_ID=$(bash "${SCRIPT_DIR}/submit_diamond_nr.sh" "${PARALLEL_ARGS[@]}" | grep -oP '(?<=job: )\d+')

FC_ARGS=( --assembly-dir "${ASM_DIR}" --annotation-dir "${ANN_DIR}" --after "${STEP2_ID}" )
[[ -n "${TEST_FLAG}" ]] && FC_ARGS+=( ${TEST_FLAG} )
STEP7_ID=$(bash "${SCRIPT_DIR}/submit_featurecounts.sh" "${FC_ARGS[@]}" | grep -oP '(?<=job: )\d+')

echo "  Step 3 KOfamScan    : ${STEP3_ID}"
echo "  Step 4 dbCAN3       : ${STEP4_ID}"
echo "  Step 5 PHI-base     : ${STEP5_ID}"
echo "  Step 6 DIAMOND NR   : ${STEP6_ID}  (slowest — ~8-16 hr/sample)"
echo "  Step 7 featureCounts: ${STEP7_ID}"
echo ""

# =============================================================================
# Summary
# =============================================================================

ALL_IDS="${STEP1_ID},${STEP2_ID},${STEP3_ID},${STEP4_ID},${STEP5_ID},${STEP6_ID},${STEP7_ID}"
[[ -n "${STEP0_ID}" ]] && ALL_IDS="${STEP0_ID},${ALL_IDS}"

echo "============================================================"
echo "All jobs submitted."
echo ""
echo "Dependency chain:"
[[ -n "${STEP0_ID}" ]] && echo "  Step 0 databases    : ${STEP0_ID}"
echo "  Step 1 Tiara        : ${STEP1_ID}"
echo "  Step 2 MetaEuk      : ${STEP2_ID}  → waits on step 1"
echo "  Steps 3-7           : parallel   → all wait on step 2"
echo ""
echo "Monitor:"
echo "  squeue -j ${ALL_IDS}"
echo "  tail -f ${ANN_DIR}/logs/metaeuk/metaeuk_${STEP2_ID}_*.out"
echo ""
echo "Logs: ${ANN_DIR}/logs/"
echo ""
echo "Step 8 — run after all jobs complete:"
echo "  Rscript ${SCRIPT_DIR}/08_integrate.R \\"
echo "    --annotation-dir ${ANN_DIR} \\"
echo "    --assembly-dir   ${ASM_DIR} \\"
echo "    --db-dir         ${DB_DIR}"
echo "============================================================"
echo ""
# =============================================================================
# RECOVERY — what to do if a step fails
# =============================================================================
# SLURM dependency type is 'afterok', meaning: if any task in an array job
# exits non-zero, all downstream dependent jobs are cancelled with status
# 'DependencyNeverSatisfied'. The idempotency checks in each job script mean
# completed samples are skipped on resubmission, so recovery is safe.
#
# To recover after a failure:
#
# 1. Find failed samples:
#      grep -l "FAILED\|Error\|error" ${ANN_DIR}/logs/<step>/*.err
#      or check which output files are missing:
#      for s in $(cat ${ANN_DIR}/sample_list.txt); do
#          [[ ! -f ${ANN_DIR}/metaeuk/${s}/${s}.fas ]] && echo "MISSING: ${s}"
#      done
#
# 2. Rerun only the failed samples (all other steps are still intact):
#      bash submit_metaeuk.sh \
#        --assembly-dir ${ASM_DIR} --annotation-dir ${ANN_DIR} \
#        --samples "FAILED_SAMPLE_01,FAILED_SAMPLE_02"
#
# 3. Once the rerun completes, resubmit the downstream steps for those samples:
#      bash submit_kofam.sh --annotation-dir ${ANN_DIR} \
#        --samples "FAILED_SAMPLE_01,FAILED_SAMPLE_02" --after <RERUN_JOB_ID>
#      # repeat for dbcan, phibase, diamond_nr, featurecounts
#
# Note: steps 3-7 are independent of each other. If only step 6 (DIAMOND NR)
# fails, steps 3-5 and 7 are unaffected and you only need to rerun step 6.
