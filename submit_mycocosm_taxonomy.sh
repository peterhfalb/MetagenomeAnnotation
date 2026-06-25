#!/bin/bash
# =============================================================================
# submit_mycocosm_taxonomy.sh
# Wrapper: submits 06c_mycocosm_taxonomy.sh as a Slurm array job.
# Runs additionally alongside submit_mmseqs_taxonomy.sh (UniRef90) — both can
# run in parallel with kofam, pfam, dbcan, phibase, and featurecounts, all
# depending only on MetaEuk (step 02) completing.
#
# Usage:
#   bash submit_mycocosm_taxonomy.sh \
#     --annotation-dir <path>
#     [--samples "S1,S2,..."]
#     [--test N]
#     [--time HH:MM:SS]
#     [--after JOB_ID]
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANNOTATION_DIR=""; TEST_N=""; SAMPLES_FILTER=""; TIME_OVERRIDE=""; AFTER_JOB=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotation-dir) ANNOTATION_DIR="$2"; shift 2 ;;
        --samples)        SAMPLES_FILTER="$2"; shift 2 ;;
        --test)           TEST_N="${2:-2}";     shift 2 ;;
        --time)           TIME_OVERRIDE="$2";  shift 2 ;;
        --after)          AFTER_JOB="$2";      shift 2 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "${ANNOTATION_DIR}" ]] && { echo "ERROR: --annotation-dir required." >&2; exit 1; }
[[ ! -d "${ANNOTATION_DIR}" ]] && { echo "ERROR: annotation dir not found: ${ANNOTATION_DIR}" >&2; exit 1; }

mkdir -p "${ANNOTATION_DIR}"/{mycocosm_taxonomy,logs/mycocosm_taxonomy}

SAMPLE_LIST="${ANNOTATION_DIR}/sample_list.txt"
[[ ! -f "${SAMPLE_LIST}" ]] && { echo "ERROR: ${SAMPLE_LIST} not found." >&2; exit 1; }

ACTIVE_LIST="${SAMPLE_LIST}"
if [[ -n "${SAMPLES_FILTER}" ]]; then
    ACTIVE_LIST="${ANNOTATION_DIR}/sample_list_mycocosm_run.txt"
    > "${ACTIVE_LIST}"
    IFS=',' read -ra NAMES <<< "${SAMPLES_FILTER}"
    for s in "${NAMES[@]}"; do
        s="${s// /}"; grep -qx "${s}" "${SAMPLE_LIST}" && echo "${s}" >> "${ACTIVE_LIST}" \
            || echo "WARNING: '${s}' not in sample list — skipping" >&2
    done
fi

N=$(wc -l < "${ACTIVE_LIST}")
[[ "${N}" -eq 0 ]] && { echo "ERROR: No valid samples." >&2; exit 1; }

if [[ -n "${TEST_N}" ]]; then
    [[ "${TEST_N}" -gt "${N}" ]] && TEST_N="${N}"
    ARRAY_RANGE="1-${TEST_N}%${TEST_N}"
    echo "TEST MODE: ${TEST_N} of ${N} samples."
else
    # CPU-bound, much smaller DB than UniRef90 — higher concurrency than the
    # UniRef90 step's GPU-node-limited %5.
    ARRAY_RANGE="1-${N}%10"
    echo "Submitting ${N} samples."
fi

SBATCH_ARGS=(
    --array="${ARRAY_RANGE}"
    --output="${ANNOTATION_DIR}/logs/mycocosm_taxonomy/myco_tax_%A_%a.out"
    --error="${ANNOTATION_DIR}/logs/mycocosm_taxonomy/myco_tax_%A_%a.err"
    --export="ANNOTATION_DIR=${ANNOTATION_DIR},ACTIVE_LIST=${ACTIVE_LIST}"
)
[[ -n "${TIME_OVERRIDE}" ]] && SBATCH_ARGS+=(--time="${TIME_OVERRIDE}")
[[ -n "${AFTER_JOB}"    ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER_JOB}")

JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" "${SCRIPT_DIR}/06c_mycocosm_taxonomy.sh" | awk '{print $NF}')
echo "Submitted Mycocosm/Phytozome taxonomy job: ${JOB_ID}"
echo "Logs: ${ANNOTATION_DIR}/logs/mycocosm_taxonomy/myco_tax_${JOB_ID}_*.out"
