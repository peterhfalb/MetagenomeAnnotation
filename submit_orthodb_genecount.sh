#!/bin/bash
# =============================================================================
# submit_orthodb_genecount.sh
# Wrapper: submits 01b_orthodb_genecount.sh as a Slurm array job.
#
# UNLIKE every other submit_*.sh in this pipeline, this step has NO
# dependency on assembly, Tiara, or MetaEuk — it only needs QC'd reads and
# the OrthoDB database (00_setup_databases.sh Section 9). It can be
# submitted immediately, in parallel with everything else, including step 1.
#
# Sample names are still resolved from ${ANNOTATION_DIR}/sample_list.txt for
# consistency with every other step (that file already exists once the
# annotation directory is initialized from the assembly output) — this
# step's own work just doesn't depend on anything else having run yet.
#
# Usage:
#   bash submit_orthodb_genecount.sh \
#     --annotation-dir <path> \
#     --qc-dir <path to QC'd filtered FASTQ directory>
#     [--samples "S1,S2,..."]
#     [--test N]
#     [--time HH:MM:SS]
#     [--after JOB_ID]   (not normally needed — this step has no real dependency)
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ANNOTATION_DIR=""; QC_DIR=""; TEST_N=""; SAMPLES_FILTER=""; TIME_OVERRIDE=""; AFTER_JOB=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotation-dir) ANNOTATION_DIR="$2"; shift 2 ;;
        --qc-dir)         QC_DIR="$2";         shift 2 ;;
        --samples)        SAMPLES_FILTER="$2"; shift 2 ;;
        --test)           TEST_N="${2:-2}";     shift 2 ;;
        --time)           TIME_OVERRIDE="$2";  shift 2 ;;
        --after)          AFTER_JOB="$2";      shift 2 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "${ANNOTATION_DIR}" ]] && { echo "ERROR: --annotation-dir required." >&2; exit 1; }
[[ -z "${QC_DIR}"         ]] && { echo "ERROR: --qc-dir required." >&2; exit 1; }
[[ ! -d "${ANNOTATION_DIR}" ]] && { echo "ERROR: annotation dir not found: ${ANNOTATION_DIR}" >&2; exit 1; }
[[ ! -d "${QC_DIR}" ]] && { echo "ERROR: QC dir not found: ${QC_DIR}" >&2; exit 1; }

mkdir -p "${ANNOTATION_DIR}"/{orthodb_genecount,logs/orthodb_genecount}

SAMPLE_LIST="${ANNOTATION_DIR}/sample_list.txt"
[[ ! -f "${SAMPLE_LIST}" ]] && { echo "ERROR: ${SAMPLE_LIST} not found." >&2; exit 1; }

ACTIVE_LIST="${SAMPLE_LIST}"
if [[ -n "${SAMPLES_FILTER}" ]]; then
    ACTIVE_LIST="${ANNOTATION_DIR}/sample_list_orthodb_run.txt"
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
    ARRAY_RANGE="1-${N}%10"
    echo "Submitting ${N} samples."
fi

SBATCH_ARGS=(
    --array="${ARRAY_RANGE}"
    --output="${ANNOTATION_DIR}/logs/orthodb_genecount/orthodb_%A_%a.out"
    --error="${ANNOTATION_DIR}/logs/orthodb_genecount/orthodb_%A_%a.err"
    --export="ANNOTATION_DIR=${ANNOTATION_DIR},QC_DIR=${QC_DIR},ACTIVE_LIST=${ACTIVE_LIST}"
)
[[ -n "${TIME_OVERRIDE}" ]] && SBATCH_ARGS+=(--time="${TIME_OVERRIDE}")
[[ -n "${AFTER_JOB}"    ]] && SBATCH_ARGS+=(--dependency="afterok:${AFTER_JOB}")

JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" "${SCRIPT_DIR}/01b_orthodb_genecount.sh" | awk '{print $NF}')
echo "Submitted OrthoDB gene count job: ${JOB_ID}"
echo "Logs: ${ANNOTATION_DIR}/logs/orthodb_genecount/orthodb_${JOB_ID}_*.out"
