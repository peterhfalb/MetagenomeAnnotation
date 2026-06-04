#!/bin/bash
# =============================================================================
# 02_metaeuk.sh
# Slurm array job: eukaryotic gene prediction with MetaEuk.
#
# MetaEuk predicts gene models on eukaryotic contigs by matching contig
# sequences against a reference protein database using MMseqs2 internally.
# It handles introns and eukaryotic gene structure that Prodigal cannot.
#
# Inputs:
#   ${ANNOTATION_DIR}/euk_contigs/${SAMPLE}_euk_contigs.fna   (from step 01)
#   ${DB_DIR}/metaeuk/fungi_refseq_db                          (pre-built MMseqs2 DB)
#
# Outputs (in ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/):
#   ${SAMPLE}.fas            predicted protein sequences    → functional annotation
#   ${SAMPLE}.codon.fas      predicted CDS nucleotide seqs → optional BLAST/taxonomy
#   ${SAMPLE}.gff            gene coordinates on contigs   → featureCounts (step 07)
#   ${SAMPLE}.headersMap.tsv maps short MetaEuk IDs to full descriptive headers
#
# Do not run directly. Submit via submit_metaeuk.sh.
# =============================================================================

#SBATCH --job-name=metaG_metaeuk
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
module load metaeuk/6-a5d39d9-gcc-8.2.0-ji6jath

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=32
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

# MetaEuk sensitivity: 7.5 is recommended for metagenome searches (range 1-7.5).
# Higher values recover more divergent/novel genes but increase runtime.
SENSITIVITY=7.5

# E-value threshold for accepting a predicted gene model.
METAEUK_EVAL=0.0001

# Minimum predicted protein length (amino acids). 20 aa = 60 nt minimum ORF.
MIN_PROTEIN_LEN=20

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

EUK_CONTIGS="${ANNOTATION_DIR}/euk_contigs/${SAMPLE}_euk_contigs.fna"
OUT_DIR="${ANNOTATION_DIR}/metaeuk/${SAMPLE}"
OUT_PREFIX="${OUT_DIR}/${SAMPLE}"
TMP_DIR="/scratch.global/${USER}/metaG_metaeuk_tmp/${SAMPLE}"

PROTEINS_FAA="${OUT_PREFIX}.fas"
GFF="${OUT_PREFIX}.gff"

echo "============================================================"
echo "Sample     : ${SAMPLE}"
echo "Job ID     : ${SLURM_JOB_ID}  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start      : $(date)"
echo "============================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "${EUK_CONTIGS}" ]]; then
    echo "ERROR: Eukaryotic contigs not found: ${EUK_CONTIGS}" >&2
    echo "  Run step 01 (submit_classify_contigs.sh) first." >&2
    exit 1
fi

# Skip if contigs file is empty (sample had no eukaryotic contigs)
if [[ ! -s "${EUK_CONTIGS}" ]]; then
    echo "WARNING: Eukaryotic contigs file is empty for ${SAMPLE} — no genes to predict."
    echo "  This sample will have no MetaEuk output. Downstream steps will skip it."
    exit 0
fi

if [[ ! -f "${DB_DIR}/metaeuk/fungi_refseq_db" ]]; then
    echo "ERROR: MetaEuk target database not found: ${DB_DIR}/metaeuk/fungi_refseq_db" >&2
    echo "  Run 00_setup_databases.sh first." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: MetaEuk easy-predict
# Runs the full MetaEuk pipeline: six-frame translation of contigs, MMseqs2
# search against the fungi reference database, and gene model construction
# accounting for introns and eukaryotic splice signals.
#
# The pre-built MMseqs2 database (fungi_refseq_db) is passed as the target;
# MetaEuk detects it is already indexed and skips DB creation.
#
# NOTE: If easy-predict fails with an error about the target database format,
# pass the FASTA directly instead:
#   Replace: ${DB_DIR}/metaeuk/fungi_refseq_db
#   With:    ${DB_DIR}/metaeuk/fungi_refseq_proteins.faa
# (MetaEuk will then build the database internally in TMP_DIR each run.)
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${PROTEINS_FAA}" && -f "${GFF}" ]]; then
    echo "[SKIP] MetaEuk outputs already exist"
else
    echo "--- Step 1: MetaEuk easy-predict ---"

    metaeuk easy-predict \
        "${EUK_CONTIGS}" \
        "${DB_DIR}/metaeuk/fungi_refseq_db" \
        "${OUT_PREFIX}" \
        "${TMP_DIR}" \
        --threads ${THREADS} \
        -s ${SENSITIVITY} \
        --metaeuk-eval ${METAEUK_EVAL} \
        --min-length ${MIN_PROTEIN_LEN}

    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary statistics
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: Gene prediction summary ---"

N_PROTEINS=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
N_EUK_CONTIGS=$(grep -c '>' "${EUK_CONTIGS}" 2>/dev/null || echo 0)
N_GFF_GENES=$(grep -v '^#' "${GFF}" 2>/dev/null | wc -l || echo 0)

echo "  Eukaryotic contigs (input) : ${N_EUK_CONTIGS}"
echo "  Predicted proteins         : ${N_PROTEINS}"
echo "  GFF feature lines          : ${N_GFF_GENES}"

# Clean up tmp directory to free scratch space
rm -rf "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key outputs:"
echo "  Proteins : ${PROTEINS_FAA}"
echo "  GFF      : ${GFF}"
echo "============================================================"
