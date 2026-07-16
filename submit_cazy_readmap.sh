#!/bin/bash
# =============================================================================
# submit_cazy_readmap.sh
# Wrapper: submits 01c_cazy_readmap.sh as a Slurm array job.
#
# Like Step 1b (OrthoDB), this step has NO dependency on assembly, Tiara, or
# MetaEuk — it only needs QC'd reads and the CAZy taxonomy database
# (00_setup_databases.sh Section 10). It can be submitted immediately, in
# parallel with everything else.
#
# Usage:
#   bash submit_cazy_readmap.sh \
#     --annotation-dir <path> \
#     --qc-dir <path to QC'd filtered FASTQ directory>
#     [--samples "S1,S2,..."]
#     [--test N]
#     [--time HH:MM:SS]
#     [--max-concurrent N]   array %N cap, default 10; use 0 for no cap
#     [--after JOB_ID]       (not normally needed — no real dependency)
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANNOTATION_DIR=""; QC_DIR=""; TEST_N=""; SAMPLES_FILTER=""; TIME_OVERRIDE=""; AFTER_JOB=""; MAX_CONCURRENT=10

while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotation-dir) ANNOTATION_DIR="$2"; shift 2 ;;
        --qc-dir)         QC_DIR="$2";         shift 2 ;;
        --samples)        SAMPLES_FILTER="$2"; shift 2 ;;
        --test)           TEST_N="${2:-2}";     shift 2 ;;
        --time)           TIME_OVERRIDE="$2";  shift 2 ;;
        --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
        --after)          AFTER_JOB="$2";      shift 2 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "${ANNOTATION_DIR}" ]] && { echo "ERROR: --annotation-dir required." >&2; exit 1; }
[[ -z "${QC_DIR}"         ]] && { echo "ERROR: --qc-dir required." >&2; exit 1; }
[[ ! -d "${ANNOTATION_DIR}" ]] && { echo "ERROR: annotation dir not found: ${ANNOTATION_DIR}" >&2; exit 1; }
[[ ! -d "${QC_DIR}" ]] && { echo "ERROR: QC dir not found: ${QC_DIR}" >&2; exit 1; }

mkdir -p "${ANNOTATION_DIR}"/{cazy_readmap,logs/cazy_readmap}

SAMPLE_LIST="${ANNOTATION_DIR}/sample_list.txt"
[[ ! -f "${SAMPLE_LIST}" ]] && { echo "ERROR: ${SAMPLE_LIST} not found." >&2; exit 1; }

ACTIVE_LIST="${SAMPLE_LIST}"
if [[ -n "${SAMPLES_FILTER}" ]]; then
    ACTIVE_LIST="${ANNOTATION_DIR}/sample_list_cazy_readmap_run.txt"
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
elif [[ "${MAX_CONCURRENT}" -eq 0 ]]; then
    ARRAY_RANGE="1-${N}"
    echo "Submitting ${N} samples, no array concurrency cap (--max-concurrent 0)."
    echo "NOTE: your account/partition's own job or CPU limits may still throttle this."
else
    ARRAY_RANGE="1-${N}%${MAX_CONCURRENT}"
    echo "Submitting ${N} samples, max ${MAX_CONCURRENT} concurrent."
fi

SBATCH_ARGS=(
    --array="${ARRAY_RANGE}"
    --output="${ANNOTATION_DIR}/logs/cazy_readmap/cazy_readmap_%A_%a.out"
    --error="${ANNOTATION_DIR}/logs/cazy_readmap/cazy_readmap_%A_%a.err"
    --export="ANNOTATION_DIR=${ANNOTATION_DIR},QC_DIR=${QC_DIR},ACTIVE_LIST=${ACTIVE_LIST}"
)
[[ -n "${TIME_OVERRIDE}" ]] && SBATCH_ARGS+=(--time="${TIME_OVERRIDE}")
[[ -n "${AFTER_JOB}"    ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER_JOB}")

JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" "${SCRIPT_DIR}/01c_cazy_readmap.sh" | awk '{print $NF}')
echo "Submitted CAZy readmap job: ${JOB_ID}"
echo "Logs: ${ANNOTATION_DIR}/logs/cazy_readmap/cazy_readmap_${JOB_ID}_*.out"
