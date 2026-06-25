#!/bin/bash
# =============================================================================
# 03b_pfam.sh
# Slurm array job: Pfam domain annotation via HMMER (hmmscan).
#
# hmmscan searches predicted proteins against the Pfam-A HMM profile database
# using each family's curated "gathering threshold" (--cut_ga) rather than a
# single fixed E-value or bitscore cutoff. This is the standard, reproducible
# approach for Pfam domain annotation — every reported hit has already passed
# that family's own threshold, so (unlike KOfamScan's detail-tsv) there is no
# sub-threshold row to filter out downstream.
#
# A protein can carry multiple distinct Pfam domains at different coordinate
# ranges (e.g. a fusion protein), so --domtblout (domain-level output) is used
# rather than --tblout (sequence-level summary only).
#
# Inputs:
#   ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas   (MetaEuk proteins)
#
# Outputs:
#   ${ANNOTATION_DIR}/pfam/${SAMPLE}_pfam_domtblout.tsv
#     Full HMMER3 domtblout: one row per domain hit, every row already passes
#     that family's gathering threshold.
#   ${ANNOTATION_DIR}/pfam/${SAMPLE}_pfam_mapper.tsv
#     Simple three-column table: gene_id, Pfam_name, Pfam_accession (version
#     suffix stripped, e.g. PF01112.21 -> PF01112). One row per domain hit —
#     a gene with multiple Pfam domains gets multiple rows.
#
# Do not run directly. Submit via submit_pfam.sh.
# =============================================================================

#SBATCH --job-name=metaG_pfam
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Conda (HMMER is in metaG_annotation environment) ──────────────────────────
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
OUT_DOMTBL="${ANNOTATION_DIR}/pfam/${SAMPLE}_pfam_domtblout.tsv"
OUT_MAPPER="${ANNOTATION_DIR}/pfam/${SAMPLE}_pfam_mapper.tsv"
OUT_LOG="${ANNOTATION_DIR}/pfam/${SAMPLE}_hmmscan.log"

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
    echo "WARNING: Protein file is empty for ${SAMPLE} — no Pfam annotations possible."
    touch "${OUT_DOMTBL}" "${OUT_MAPPER}"
    exit 0
fi
if [[ ! -f "${DB_DIR}/pfam/Pfam-A.hmm.h3i" ]]; then
    echo "ERROR: Pfam-A database not found in ${DB_DIR}/pfam/" >&2
    echo "  Run 00_setup_databases.sh first." >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_DOMTBL}")"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: hmmscan — query=protein, target=Pfam-A.hmm, per-family gathering
# threshold (--cut_ga) instead of a fixed E-value. --noali suppresses full
# alignment dumps to keep output size manageable (same rationale as why
# KOfam/dbCAN outputs avoid verbose alignment text).
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_DOMTBL}" && -s "${OUT_DOMTBL}" ]]; then
    echo "[SKIP] hmmscan output already exists"
else
    echo "--- Step 1: hmmscan vs Pfam-A (--cut_ga) ---"
    hmmscan \
        --cut_ga \
        --cpu ${THREADS} \
        --domtblout "${OUT_DOMTBL}" \
        --noali \
        "${DB_DIR}/pfam/Pfam-A.hmm" \
        "${PROTEINS_FAA}" \
        > "${OUT_LOG}"
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Derive a simple gene_id -> Pfam mapper table from the domtblout.
# domtblout is a fixed-width, whitespace-delimited, '#'-commented format:
#   col 1  = target name   (Pfam family name, e.g. "Asparaginase")
#   col 2  = target acc.   (Pfam accession + version, e.g. "PF01112.21")
#   col 4  = query name    (gene_id)
# Every row already passed --cut_ga, so unlike KOfam's detail-tsv there is no
# significance flag to filter on here.
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: Extracting gene_id -> Pfam mapper table ---"

awk '!/^#/ {print $4 "\t" $1 "\t" $2}' "${OUT_DOMTBL}" \
    | sed -E 's/\t(PF[0-9]+)\.[0-9]+$/\t\1/' \
    > "${OUT_MAPPER}"

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
N_HITS=$(wc -l < "${OUT_MAPPER}")
echo "  Proteins searched   : ${N_TOTAL}"
echo "  Pfam domain hits     : ${N_HITS} (gene-domain pairs, all passing --cut_ga)"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key outputs:"
echo "  Full domtblout : ${OUT_DOMTBL}"
echo "  Mapper table   : ${OUT_MAPPER}"
echo "============================================================"
