#!/bin/bash
# =============================================================================
# 01_classify_contigs.sh
# Slurm array job: classify contigs by domain and extract eukaryotic subset.
#
# Steps:
#   1. Tiara  — deep-learning contig classifier; labels each contig as
#               eukaryota / prokarya / mitochondria / plastid / unknown
#   2. Extract — split contigs.fna into euk_contigs.fna and prok_contigs.fna
#                using samtools faidx for downstream gene prediction
#
# Inputs (from MetagenomeAssembly output):
#   ${ASSEMBLY_DIR}/assemblies/contigs/${SAMPLE}_contigs.fna
#
# Outputs:
#   ${ANNOTATION_DIR}/tiara/${SAMPLE}_tiara.txt        full Tiara classification table
#   ${ANNOTATION_DIR}/tiara/${SAMPLE}_euk_list.txt     eukaryotic contig names
#   ${ANNOTATION_DIR}/tiara/${SAMPLE}_prok_list.txt    prokaryotic contig names
#   ${ANNOTATION_DIR}/euk_contigs/${SAMPLE}_euk_contigs.fna
#   ${ANNOTATION_DIR}/prok_contigs/${SAMPLE}_prok_contigs.fna
#   ${ANNOTATION_DIR}/tiara/${SAMPLE}_tiara_summary.txt  per-class contig counts
#
# Do not run directly. Submit via submit_classify_contigs.sh.
# =============================================================================

#SBATCH --job-name=metaG_tiara
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Conda (Tiara is in metaG_tiara environment) ───────────────────────────────
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_tiara

# ── Modules ───────────────────────────────────────────────────────────────────
module load samtools/1.21

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=16

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

CONTIGS_FNA="${ASSEMBLY_DIR}/assemblies/contigs/${SAMPLE}_contigs.fna"

TIARA_OUT="${ANNOTATION_DIR}/tiara/${SAMPLE}_tiara.txt"
EUK_LIST="${ANNOTATION_DIR}/tiara/${SAMPLE}_euk_list.txt"
PROK_LIST="${ANNOTATION_DIR}/tiara/${SAMPLE}_prok_list.txt"
EUK_CONTIGS="${ANNOTATION_DIR}/euk_contigs/${SAMPLE}_euk_contigs.fna"
PROK_CONTIGS="${ANNOTATION_DIR}/prok_contigs/${SAMPLE}_prok_contigs.fna"
SUMMARY="${ANNOTATION_DIR}/tiara/${SAMPLE}_tiara_summary.txt"

echo "============================================================"
echo "Sample     : ${SAMPLE}"
echo "Job ID     : ${SLURM_JOB_ID}  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start      : $(date)"
echo "============================================================"

# ── Validate input ────────────────────────────────────────────────────────────
if [[ ! -f "${CONTIGS_FNA}" ]]; then
    echo "ERROR: Contigs file not found: ${CONTIGS_FNA}" >&2
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Tiara contig classification
# Classifies each contig by k-mer composition using a deep bidirectional LSTM.
# first_classifier: eukaryota | prokarya | mitochondria | plastid | unknown
# second_classifier: for prokarya only → bacteria | archaea
# --probabilities: adds per-class probability columns (kept for QC)
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${TIARA_OUT}" ]]; then
    echo "[SKIP] Tiara classification already exists"
else
    echo "--- Step 1: Tiara classification ---"
    tiara \
        -i "${CONTIGS_FNA}" \
        -o "${TIARA_OUT}" \
        --probabilities \
        -t ${THREADS} \
        --min_len 500
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Extract contig name lists by domain
# Tiara output columns: sequence_name, first_classifier, second_classifier,
#   [per-class probabilities if --probabilities used]
# Mitochondria and plastid contigs are excluded from both lists — they have
# very different gene structure and would confuse both MetaEuk and Prodigal.
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: Extracting contig name lists ---"

awk '$2 == "eukaryota"  {print $1}' "${TIARA_OUT}" > "${EUK_LIST}"
awk '$2 == "prokarya"   {print $1}' "${TIARA_OUT}" > "${PROK_LIST}"

# Summary counts (written to file and stdout for logs)
N_TOTAL=$(grep -c '>' "${CONTIGS_FNA}")
N_EUK=$(wc -l < "${EUK_LIST}")
N_PROK=$(wc -l < "${PROK_LIST}")
N_MITO=$(awk '$2 == "mitochondria" {c++} END {print c+0}' "${TIARA_OUT}")
N_PLASTID=$(awk '$2 == "plastid"   {c++} END {print c+0}' "${TIARA_OUT}")
N_UNKNOWN=$(awk '$2 == "unknown"   {c++} END {print c+0}' "${TIARA_OUT}")

pct() { awk "BEGIN{printf \"%.1f\", ${N_TOTAL}>0 ? ${1}/${N_TOTAL}*100 : 0}"; }
{
    echo "Sample          : ${SAMPLE}"
    echo "Total contigs   : ${N_TOTAL}"
    echo "  Eukaryota     : ${N_EUK}  ($(pct "${N_EUK}")%)"
    echo "  Prokarya      : ${N_PROK}  ($(pct "${N_PROK}")%)"
    echo "  Mitochondria  : ${N_MITO}"
    echo "  Plastid       : ${N_PLASTID}"
    echo "  Unknown       : ${N_UNKNOWN}"
} | tee "${SUMMARY}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Extract eukaryotic and prokaryotic contigs as separate FASTAs
# samtools faidx -r reads sequence names from a file, avoids shell arg limits.
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 3: Extracting FASTA subsets ---"

# Create faidx index if not already present
[[ ! -f "${CONTIGS_FNA}.fai" ]] && samtools faidx "${CONTIGS_FNA}"

if [[ -f "${EUK_CONTIGS}" ]]; then
    echo "[SKIP] Eukaryotic contigs FASTA already exists"
else
    if [[ "${N_EUK}" -gt 0 ]]; then
        samtools faidx "${CONTIGS_FNA}" \
            -r "${EUK_LIST}" \
            -o "${EUK_CONTIGS}"
        echo "  Extracted ${N_EUK} eukaryotic contigs"
    else
        echo "WARNING: No eukaryotic contigs found for ${SAMPLE} — creating empty file"
        touch "${EUK_CONTIGS}"
    fi
fi

if [[ -f "${PROK_CONTIGS}" ]]; then
    echo "[SKIP] Prokaryotic contigs FASTA already exists"
else
    if [[ "${N_PROK}" -gt 0 ]]; then
        samtools faidx "${CONTIGS_FNA}" \
            -r "${PROK_LIST}" \
            -o "${PROK_CONTIGS}"
        echo "  Extracted ${N_PROK} prokaryotic contigs"
    else
        echo "WARNING: No prokaryotic contigs found for ${SAMPLE} — creating empty file"
        touch "${PROK_CONTIGS}"
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "============================================================"
