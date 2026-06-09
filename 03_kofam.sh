#!/bin/bash
# =============================================================================
# 03_kofam.sh
# Slurm array job: KEGG KO assignment via KOfamScan.
#
# KOfamScan searches predicted proteins against KEGG KOfam HMM profiles
# (one profile per KO number). Each hit above the KO-specific score threshold
# is assigned that KO. No KEGG license required — KOfam profiles are freely
# distributed by GenomeNet.
#
# Inputs:
#   ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas   (MetaEuk proteins)
#
# Outputs:
#   ${ANNOTATION_DIR}/kofam/${SAMPLE}_kofam.tsv
#     detail-tsv format: tab-separated, one row per hit (significant and not).
#     Column 1: '*' if hit is above threshold (significant), else blank.
#     Filter to '*' rows in R for downstream analysis.
#   ${ANNOTATION_DIR}/kofam/${SAMPLE}_kofam_mapper.tsv
#     Simple two-column table: gene_id → KO (significant hits only, one per gene).
#     Used as a quick lookup table for pathway mapping.
#
# Do not run directly. Submit via submit_kofam.sh.
# =============================================================================

#SBATCH --job-name=metaG_kofam
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Conda (KOfamScan is in metaG_annotation environment) ─────────────────────
set +u
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_annotation
set -u

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=32
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

PROTEINS_FAA="${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas"
OUT_TSV="${ANNOTATION_DIR}/kofam/${SAMPLE}_kofam.tsv"
OUT_MAPPER="${ANNOTATION_DIR}/kofam/${SAMPLE}_kofam_mapper.tsv"
TMP_DIR="/scratch.global/${USER}/metaG_kofam_tmp/${SAMPLE}"

echo "============================================================"
echo "Sample     : ${SAMPLE}"
echo "Job ID     : ${SLURM_JOB_ID}  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start      : $(date)"
echo "============================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "${PROTEINS_FAA}" ]]; then
    echo "ERROR: MetaEuk proteins not found: ${PROTEINS_FAA}" >&2
    echo "  Run step 02 (submit_metaeuk.sh) first." >&2
    exit 1
fi
if [[ ! -s "${PROTEINS_FAA}" ]]; then
    echo "WARNING: Protein file is empty for ${SAMPLE} — no KO annotations possible."
    touch "${OUT_TSV}" "${OUT_MAPPER}"
    exit 0
fi
if [[ ! -d "${DB_DIR}/kofam/profiles" || ! -f "${DB_DIR}/kofam/ko_list" ]]; then
    echo "ERROR: KOfam database not found in ${DB_DIR}/kofam/" >&2
    echo "  Run 00_setup_databases.sh first." >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_TSV}")" "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: KOfamScan — annotate proteins with KEGG KO numbers
# detail-tsv format captures all hits (above and below threshold) with scores.
# Column 1 == '*' marks hits above the KO-specific threshold (significant).
# Keeping sub-threshold hits allows re-filtering with different cutoffs in R.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_TSV}" && -s "${OUT_TSV}" ]]; then
    echo "[SKIP] KOfamScan output already exists"
else
    echo "--- Step 1: KOfamScan ---"
    exec_annotation \
        -f detail-tsv \
        -o "${OUT_TSV}" \
        --ko-list "${DB_DIR}/kofam/ko_list" \
        --profile "${DB_DIR}/kofam/profiles" \
        --cpu ${THREADS} \
        --tmp-dir "${TMP_DIR}" \
        "${PROTEINS_FAA}"
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Extract significant hits as a simple gene→KO mapper table
# Filters to '*' (threshold-passing) rows, outputs gene_id and KO only.
# When a gene has multiple significant KO hits, all are kept (one per line).
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: Extracting significant hits ---"

awk 'NF >= 3 && $1 == "*" {print $2 "\t" $3}' "${OUT_TSV}" > "${OUT_MAPPER}"

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
N_SIG=$(wc -l < "${OUT_MAPPER}")
echo "  Proteins searched      : ${N_TOTAL}"
echo "  Significant KO hits    : ${N_SIG} (gene-KO pairs above threshold)"

rm -rf "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key outputs:"
echo "  Full results : ${OUT_TSV}"
echo "  Mapper table : ${OUT_MAPPER}"
echo "============================================================"
