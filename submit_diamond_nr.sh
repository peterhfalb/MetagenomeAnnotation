#!/bin/bash
# =============================================================================
# submit_diamond_nr.sh
# Wrapper: submits 06_diamond_nr.sh as a Slurm array job.
# This is the slowest annotation step (~8-16 hr per sample).
# Can run in parallel with kofam, dbcan, phibase after MetaEuk completes.
# Max concurrency is kept lower (3 jobs) to avoid NR database I/O contention.
#
# Usage:
#   bash submit_diamond_nr.sh \
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

mkdir -p "${ANNOTATION_DIR}"/{diamond_nr,logs/diamond_nr}

SAMPLE_LIST="${ANNOTATION_DIR}/sample_list.txt"
[[ ! -f "${SAMPLE_LIST}" ]] && { echo "ERROR: ${SAMPLE_LIST} not found." >&2; exit 1; }

ACTIVE_LIST="${SAMPLE_LIST}"
if [[ -n "${SAMPLES_FILTER}" ]]; then
    ACTIVE_LIST="${ANNOTATION_DIR}/sample_list_diamond_run.txt"
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
    # Low concurrency: NR database reads can saturate shared storage at high concurrency
    ARRAY_RANGE="1-${N}%3"
    echo "Submitting ${N} samples (max 3 concurrent — NR I/O throttle)."
fi

SBATCH_ARGS=(
    --array="${ARRAY_RANGE}"
    --output="${ANNOTATION_DIR}/logs/diamond_nr/diamond_nr_%A_%a.out"
    --error="${ANNOTATION_DIR}/logs/diamond_nr/diamond_nr_%A_%a.err"
    --export="ANNOTATION_DIR=${ANNOTATION_DIR},ACTIVE_LIST=${ACTIVE_LIST}"
)
[[ -n "${TIME_OVERRIDE}" ]] && SBATCH_ARGS+=(--time="${TIME_OVERRIDE}")
[[ -n "${AFTER_JOB}"    ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER_JOB}")

JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" "${SCRIPT_DIR}/06_diamond_nr.sh" | awk '{print $NF}')
echo "Submitted DIAMOND NR job: ${JOB_ID}"
echo "Logs: ${ANNOTATION_DIR}/logs/diamond_nr/diamond_nr_${JOB_ID}_*.out"
echo ""
echo "NOTE: Expect ~8-16 hr per sample. Monitor with: squeue -j ${JOB_ID}"
