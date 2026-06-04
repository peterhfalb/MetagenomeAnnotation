#!/bin/bash
# =============================================================================
# submit_metaeuk.sh
# Wrapper: submits 02_metaeuk.sh as a Slurm array job.
# Run from the login node after step 01 (Tiara) is complete.
#
# Usage:
#   bash submit_metaeuk.sh \
#     --assembly-dir  <path>    same assembly output dir used in step 01
#     --annotation-dir <path>   same annotation dir used in step 01
#     [--samples "S1,S2,..."    rerun specific samples only]
#     [--test N                 run only the first N samples]
#     [--time HH:MM:SS          override default walltime]
#     [--after JOB_ID           hold until Tiara job JOB_ID completes]
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Parse arguments ───────────────────────────────────────────────────────────
ASSEMBLY_DIR=""
ANNOTATION_DIR=""
TEST_N=""
SAMPLES_FILTER=""
TIME_OVERRIDE=""
AFTER_JOB=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly-dir)   ASSEMBLY_DIR="$2";    shift 2 ;;
        --annotation-dir) ANNOTATION_DIR="$2";  shift 2 ;;
        --samples)        SAMPLES_FILTER="$2";  shift 2 ;;
        --test)           TEST_N="${2:-2}";      shift 2 ;;
        --time)           TIME_OVERRIDE="$2";   shift 2 ;;
        --after)          AFTER_JOB="$2";       shift 2 ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            echo "Usage: bash submit_metaeuk.sh --assembly-dir <path> --annotation-dir <path> [--samples \"S1,S2\"] [--test N] [--time HH:MM:SS] [--after JOB_ID]" >&2
            exit 1 ;;
    esac
done

if [[ -z "${ASSEMBLY_DIR}" || -z "${ANNOTATION_DIR}" ]]; then
    echo "ERROR: --assembly-dir and --annotation-dir are both required." >&2
    exit 1
fi
if [[ ! -d "${ANNOTATION_DIR}" ]]; then
    echo "ERROR: --annotation-dir does not exist: ${ANNOTATION_DIR}" >&2
    echo "  Run submit_classify_contigs.sh first." >&2
    exit 1
fi

# ── Create output directories ─────────────────────────────────────────────────
mkdir -p "${ANNOTATION_DIR}"/{metaeuk,logs/metaeuk}

# ── Resolve sample list ───────────────────────────────────────────────────────
SAMPLE_LIST="${ANNOTATION_DIR}/sample_list.txt"
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found at ${SAMPLE_LIST}" >&2
    echo "  Expected to be created by submit_classify_contigs.sh." >&2
    exit 1
fi

TOTAL=$(wc -l < "${SAMPLE_LIST}")
ACTIVE_LIST="${SAMPLE_LIST}"

if [[ -n "${SAMPLES_FILTER}" ]]; then
    ACTIVE_LIST="${ANNOTATION_DIR}/sample_list_run.txt"
    > "${ACTIVE_LIST}"
    IFS=',' read -ra SAMPLE_NAMES <<< "${SAMPLES_FILTER}"
    for s in "${SAMPLE_NAMES[@]}"; do
        s="${s// /}"
        if grep -qx "${s}" "${SAMPLE_LIST}"; then
            echo "${s}" >> "${ACTIVE_LIST}"
        else
            echo "WARNING: '${s}' not found in sample list — skipping" >&2
        fi
    done
fi

N=$(wc -l < "${ACTIVE_LIST}")
if [[ "${N}" -eq 0 ]]; then
    echo "ERROR: No valid samples to submit." >&2; exit 1
fi

# ── Determine array range ─────────────────────────────────────────────────────
if [[ -n "${TEST_N}" ]]; then
    [[ "${TEST_N}" -gt "${N}" ]] && TEST_N="${N}"
    ARRAY_RANGE="1-${TEST_N}%${TEST_N}"
    echo "TEST MODE: submitting first ${TEST_N} of ${N} samples."
else
    # Max 5 concurrent — MetaEuk is memory-heavy; avoid saturating the cluster
    ARRAY_RANGE="1-${N}%5"
    echo "Submitting ${N} sample(s) (${TOTAL} total)."
fi

echo "Active sample list : ${ACTIVE_LIST}"
echo "Array range        : ${ARRAY_RANGE}"
[[ -n "${TIME_OVERRIDE}" ]] && echo "Walltime override  : ${TIME_OVERRIDE}"
[[ -n "${AFTER_JOB}"    ]] && echo "Dependency         : afterok:${AFTER_JOB}"
echo "Submitting..."

# ── Submit ────────────────────────────────────────────────────────────────────
SBATCH_ARGS=(
    --array="${ARRAY_RANGE}"
    --output="${ANNOTATION_DIR}/logs/metaeuk/metaeuk_%A_%a.out"
    --error="${ANNOTATION_DIR}/logs/metaeuk/metaeuk_%A_%a.err"
    --export="ASSEMBLY_DIR=${ASSEMBLY_DIR},ANNOTATION_DIR=${ANNOTATION_DIR},ACTIVE_LIST=${ACTIVE_LIST}"
)

[[ -n "${TIME_OVERRIDE}" ]] && SBATCH_ARGS+=(--time="${TIME_OVERRIDE}")
[[ -n "${AFTER_JOB}"    ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER_JOB}")

JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" "${SCRIPT_DIR}/02_metaeuk.sh" | awk '{print $NF}')

echo "Submitted MetaEuk array job: ${JOB_ID}"
echo "Logs: ${ANNOTATION_DIR}/logs/metaeuk/metaeuk_${JOB_ID}_*.out"
echo ""
echo "When complete, submit functional annotation steps:"
echo "  bash submit_kofam.sh    --annotation-dir ${ANNOTATION_DIR} --after ${JOB_ID}"
echo "  bash submit_dbcan.sh    --annotation-dir ${ANNOTATION_DIR} --after ${JOB_ID}"
echo "  bash submit_phibase.sh  --annotation-dir ${ANNOTATION_DIR} --after ${JOB_ID}"
echo "  bash submit_diamond_nr.sh --annotation-dir ${ANNOTATION_DIR} --after ${JOB_ID}"
