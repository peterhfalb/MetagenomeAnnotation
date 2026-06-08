#!/bin/bash
# =============================================================================
# submit_mmseqs_taxonomy.sh
# Wrapper: submits 06_mmseqs_taxonomy.sh as a Slurm array job.
# Requests GPU nodes (a100 partition). Can run in parallel with kofam, dbcan,
# phibase, and featurecounts after MetaEuk (step 02) completes.
#
# Usage:
#   bash submit_mmseqs_taxonomy.sh \
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

mkdir -p "${ANNOTATION_DIR}"/{mmseqs_taxonomy,logs/mmseqs_taxonomy}

SAMPLE_LIST="${ANNOTATION_DIR}/sample_list.txt"
[[ ! -f "${SAMPLE_LIST}" ]] && { echo "ERROR: ${SAMPLE_LIST} not found." >&2; exit 1; }

ACTIVE_LIST="${SAMPLE_LIST}"
if [[ -n "${SAMPLES_FILTER}" ]]; then
    ACTIVE_LIST="${ANNOTATION_DIR}/sample_list_mmseqs_run.txt"
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
    # GPU jobs: run up to 5 concurrently (GPU nodes are limited)
    ARRAY_RANGE="1-${N}%5"
    echo "Submitting ${N} samples (max 5 concurrent — GPU node limit)."
fi

SBATCH_ARGS=(
    --array="${ARRAY_RANGE}"
    --output="${ANNOTATION_DIR}/logs/mmseqs_taxonomy/mmseqs_tax_%A_%a.out"
    --error="${ANNOTATION_DIR}/logs/mmseqs_taxonomy/mmseqs_tax_%A_%a.err"
    --export="ANNOTATION_DIR=${ANNOTATION_DIR},ACTIVE_LIST=${ACTIVE_LIST}"
)
[[ -n "${TIME_OVERRIDE}" ]] && SBATCH_ARGS+=(--time="${TIME_OVERRIDE}")
[[ -n "${AFTER_JOB}"    ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER_JOB}")

JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" "${SCRIPT_DIR}/06_mmseqs_taxonomy.sh" | awk '{print $NF}')
echo "Submitted MMseqs2 taxonomy job: ${JOB_ID}"
echo "Logs: ${ANNOTATION_DIR}/logs/mmseqs_taxonomy/mmseqs_tax_${JOB_ID}_*.out"
