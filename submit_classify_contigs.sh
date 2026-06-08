#!/bin/bash
# =============================================================================
# submit_classify_contigs.sh
# Wrapper: creates output directories, resolves sample list, submits array job.
# Run from the login node after database setup is complete.
#
# Usage:
#   bash submit_classify_contigs.sh \
#     --assembly-dir  <path>    output directory from MetagenomeAssembly pipeline
#     --annotation-dir <path>   annotation output directory (will be created)
#     [--samples "S1,S2,..."    run only these named samples]
#     [--test N                 run only the first N samples (default 2)]
#     [--time HH:MM:SS          override default walltime]
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
            echo "Usage: bash submit_classify_contigs.sh --assembly-dir <path> --annotation-dir <path> [--samples \"S1,S2\"] [--test N] [--time HH:MM:SS] [--after JOB_ID]" >&2
            exit 1 ;;
    esac
done

if [[ -z "${ASSEMBLY_DIR}" || -z "${ANNOTATION_DIR}" ]]; then
    echo "ERROR: --assembly-dir and --annotation-dir are both required." >&2
    exit 1
fi
if [[ ! -d "${ASSEMBLY_DIR}" ]]; then
    echo "ERROR: --assembly-dir does not exist: ${ASSEMBLY_DIR}" >&2
    exit 1
fi

# ── Create output directory structure ────────────────────────────────────────
mkdir -p "${ANNOTATION_DIR}"/{tiara,euk_contigs,prok_contigs,logs/tiara}

# ── Build sample list from assembly output ────────────────────────────────────
ASSEMBLY_SAMPLE_LIST="${ASSEMBLY_DIR}/sample_list.txt"
if [[ ! -f "${ASSEMBLY_SAMPLE_LIST}" ]]; then
    echo "ERROR: Sample list not found at ${ASSEMBLY_SAMPLE_LIST}" >&2
    echo "  Expected: the sample_list.txt created by submit_assembly.sh" >&2
    exit 1
fi

cp "${ASSEMBLY_SAMPLE_LIST}" "${ANNOTATION_DIR}/sample_list.txt"
SAMPLE_LIST="${ANNOTATION_DIR}/sample_list.txt"
TOTAL=$(wc -l < "${SAMPLE_LIST}")

# ── Filter to specific samples if --samples provided ─────────────────────────
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
    ARRAY_RANGE="1-${N}%10"
    echo "Submitting ${N} sample(s) (${TOTAL} total in assembly dir)."
fi

echo "Active sample list : ${ACTIVE_LIST}"
echo "Array range        : ${ARRAY_RANGE}"
[[ -n "${TIME_OVERRIDE}" ]] && echo "Walltime override  : ${TIME_OVERRIDE}"
echo "Submitting..."

# ── Submit ────────────────────────────────────────────────────────────────────
SBATCH_ARGS=(
    --array="${ARRAY_RANGE}"
    --output="${ANNOTATION_DIR}/logs/tiara/tiara_%A_%a.out"
    --error="${ANNOTATION_DIR}/logs/tiara/tiara_%A_%a.err"
    --export="ASSEMBLY_DIR=${ASSEMBLY_DIR},ANNOTATION_DIR=${ANNOTATION_DIR},ACTIVE_LIST=${ACTIVE_LIST}"
)

[[ -n "${TIME_OVERRIDE}" ]] && SBATCH_ARGS+=(--time="${TIME_OVERRIDE}")
[[ -n "${AFTER_JOB}"    ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER_JOB}")

JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" "${SCRIPT_DIR}/01_classify_contigs.sh" | awk '{print $NF}')

echo "Submitted Tiara array job: ${JOB_ID}"
echo "Logs: ${ANNOTATION_DIR}/logs/tiara/tiara_${JOB_ID}_*.out"
echo ""
echo "When complete, submit MetaEuk gene prediction with:"
echo "  bash submit_metaeuk.sh --assembly-dir ${ASSEMBLY_DIR} --annotation-dir ${ANNOTATION_DIR}"
