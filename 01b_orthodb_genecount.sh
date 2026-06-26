#!/bin/bash
# =============================================================================
# 01b_orthodb_genecount.sh
# Slurm array job: genome-equivalent estimation via DIAMOND vs OrthoDB Fungi
# near-single-copy orthologs.
#
# Provides an ADDITIVE alternative to the Asparaginase (Pfam PF01112)
# pseudo-genome normalization already computed in 08_integrate.R: instead of
# one marker gene, this maps raw QC'd reads directly against near-single-copy
# Fungi-level orthologs from OrthoDB v12 (NCBI taxid 4751 — broader than
# Dikarya, which OrthoDB doesn't compute orthologous groups for at all; see
# 00_setup_databases.sh Section 9). 08_integrate.R later averages read counts
# per ortholog group (to avoid double-counting from within-group sequence
# redundancy) and takes the geometric mean across all groups — far more
# robust than relying on a single marker gene, since any one ortholog's
# noise/absence is washed out by averaging over many.
#
# UNLIKE every other step in this pipeline, this one operates on raw QC'd
# reads, PRE-ASSEMBLY — it does not need Tiara, MetaEuk, or any assembled
# contigs, and has NO dependency on any other step. It can run immediately,
# as soon as the OrthoDB database (00_setup_databases.sh Section 9) exists.
#
# Inputs:
#   ${QC_DIR}/${SAMPLE}_R1_filtered.fastq.gz
#   ${QC_DIR}/${SAMPLE}_R2_filtered.fastq.gz
#
# Outputs:
#   ${ANNOTATION_DIR}/orthodb_genecount/${SAMPLE}_hits.tsv
#     Merged DIAMOND blastx hits (R1 + R2) vs the OrthoDB Fungi database.
#     Columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart,
#     qend, sstart, send, evalue, bitscore. sseqid is the OrthoDB sequence
#     ID — 08_integrate.R maps this to an ortholog group via gene2og.tsv
#     (built in 00_setup_databases.sh Section 9).
#
# Do not run directly. Submit via submit_orthodb_genecount.sh.
# =============================================================================

#SBATCH --job-name=metaG_orthodb_genecount
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=20:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

# NOTE on resources: observed runtime on one real test sample — R1 alone took
# ~4h23m (16 CPU, --sensitive mode), so R1+R2 combined is realistically
# ~8-9 hours, not the "much faster than 12 hours" originally guessed here.
# Bumped to 20 hours for safety margin across samples with larger/more
# complex read sets; tighten this once more samples have actually run.

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
module load diamond/2.0.15-gcc-8.2.0-gkldzx7

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=16
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"
ORTHODB_DB="${DB_DIR}/orthodb/fungi_orthologs"

# E-value cutoff and sensitivity matching this pipeline's existing DIAMOND
# convention (05_phibase.sh). --max-target-seqs 1: best hit per read only,
# so each read contributes to at most one OrthoDB sequence — avoids
# read-level multi-counting before the OG-level averaging step in R.
EVALUE=1e-5

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

R1_FASTQ="${QC_DIR}/${SAMPLE}_R1_filtered.fastq.gz"
R2_FASTQ="${QC_DIR}/${SAMPLE}_R2_filtered.fastq.gz"
OUT_DIR="${ANNOTATION_DIR}/orthodb_genecount"
OUT_HITS="${OUT_DIR}/${SAMPLE}_hits.tsv"
OUT_R1="${OUT_DIR}/${SAMPLE}_R1_hits.tsv"
OUT_R2="${OUT_DIR}/${SAMPLE}_R2_hits.tsv"

echo "============================================================"
echo "Sample     : ${SAMPLE}"
echo "Job ID     : ${SLURM_JOB_ID}  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start      : $(date)"
echo "============================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "${R1_FASTQ}" || ! -f "${R2_FASTQ}" ]]; then
    echo "ERROR: QC'd reads not found:" >&2
    echo "  ${R1_FASTQ}" >&2
    echo "  ${R2_FASTQ}" >&2
    exit 1
fi
if [[ ! -f "${ORTHODB_DB}.dmnd" ]]; then
    echo "ERROR: OrthoDB Fungi database not found: ${ORTHODB_DB}.dmnd" >&2
    echo "  Run 00_setup_databases.sh first (Section 9)." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: DIAMOND blastx vs OrthoDB Fungi orthologs, R1 and R2 separately
# (matching the source paper's "forward and reverse reads were mapped"),
# then merge. blastx (not blastp) because the input is raw nucleotide reads,
# not predicted proteins.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_HITS}" && -s "${OUT_HITS}" ]]; then
    echo "[SKIP] OrthoDB gene count output already exists"
else
    # Per-direction idempotency: R1 alone takes ~4+ hours, so if the job is
    # killed mid-R2 (timeout, node failure) and resubmitted, don't throw away
    # R1's already-completed work — only the top-level OUT_HITS check existed
    # before, which would have silently redone R1 from scratch on resubmit.
    if [[ -f "${OUT_R1}" && -s "${OUT_R1}" ]]; then
        echo "--- Step 1a: [SKIP] R1 output already exists ---"
    else
        echo "--- Step 1a: DIAMOND blastx vs OrthoDB (R1) ---"
        diamond blastx \
            --query "${R1_FASTQ}" \
            --db "${ORTHODB_DB}" \
            --out "${OUT_R1}" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --evalue ${EVALUE} \
            --max-target-seqs 1 \
            --sensitive \
            --threads ${THREADS}
        echo "Step 1a done: $(date)"
    fi

    if [[ -f "${OUT_R2}" && -s "${OUT_R2}" ]]; then
        echo "--- Step 1b: [SKIP] R2 output already exists ---"
    else
        echo "--- Step 1b: DIAMOND blastx vs OrthoDB (R2) ---"
        diamond blastx \
            --query "${R2_FASTQ}" \
            --db "${ORTHODB_DB}" \
            --out "${OUT_R2}" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --evalue ${EVALUE} \
            --max-target-seqs 1 \
            --sensitive \
            --threads ${THREADS}
        echo "Step 1b done: $(date)"
    fi

    cat "${OUT_R1}" "${OUT_R2}" > "${OUT_HITS}"
    rm -f "${OUT_R1}" "${OUT_R2}"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# ─────────────────────────────────────────────────────────────────────────────

N_HITS=$(wc -l < "${OUT_HITS}")
N_UNIQUE_OGSEQS=$(cut -f2 "${OUT_HITS}" | sort -u | wc -l)

echo "  Total read hits (R1+R2)        : ${N_HITS}"
echo "  Unique OrthoDB sequences hit    : ${N_UNIQUE_OGSEQS}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OUT_HITS}"
echo "============================================================"
