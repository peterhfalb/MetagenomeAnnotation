#!/bin/bash
# =============================================================================
# 01c_cazy_readmap.sh
# Slurm array job: direct read mapping against CAZy database (all kingdoms)
# with NCBI taxonomy, following Bahram 2018 (Nature).
#
# Maps raw QC'd reads (pre-assembly) against a taxonomy-enabled CAZy DIAMOND
# database (built in 00_setup_databases.sh Section 10). Unlike the assembly-
# based CAZyme annotation (MetaEuk → dbCAN), this captures reads that never
# assembled or landed on contigs too short for Tiara/MetaEuk. Maps competitively
# against all kingdoms — kingdom assignment (Fungi / Bacteria / etc.) is done
# in 08_integrate.R using the staxids column.
#
# Inputs:
#   ${QC_DIR}/${SAMPLE}_R1_filtered.fastq.gz
#   ${QC_DIR}/${SAMPLE}_R2_filtered.fastq.gz
#
# Outputs:
#   ${ANNOTATION_DIR}/cazy_readmap/${SAMPLE}_hits.tsv
#     Merged DIAMOND blastx hits (R1 + R2) vs CAZy_taxonomy database.
#     Columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart,
#     qend, sstart, send, evalue, bitscore, staxids
#     sseqid format: ACCESSION|CAZyFamily (e.g. WP_123456.1|GH5)
#     staxids: NCBI taxid of best-hit sequence — used for kingdom assignment
#
# Note: raw hits use e-value ≤ 1e-5 (DIAMOND filter). The final identity
# (≥50%) and e-value (≤ 1e-9) filters from Bahram 2018 are applied in
# 08_integrate.R at integration time, preserving the raw hits for sensitivity
# analysis with different thresholds.
#
# Do not run directly. Submit via submit_cazy_readmap.sh.
# =============================================================================

#SBATCH --job-name=metaG_cazy_readmap
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

# NOTE on time: CAZy_taxonomy.dmnd has 4M sequences — much larger than OrthoDB.
# Empirically, R1 alone takes ~3-4 hours per sample; R2 similar. Samples with
# larger read files (up to 7 GB) can approach 12h total. 24h gives safe buffer.

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
module load diamond/2.0.15-gcc-8.2.0-gkldzx7

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=16
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"
CAZY_DB="${DB_DIR}/cazy_readmap/CAZy_taxonomy"
EVALUE=1e-5

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

R1_FASTQ="${QC_DIR}/${SAMPLE}_R1_filtered.fastq.gz"
R2_FASTQ="${QC_DIR}/${SAMPLE}_R2_filtered.fastq.gz"
OUT_DIR="${ANNOTATION_DIR}/cazy_readmap"
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
if [[ ! -f "${CAZY_DB}.dmnd" ]]; then
    echo "ERROR: CAZy taxonomy database not found: ${CAZY_DB}.dmnd" >&2
    echo "  Run 00_setup_databases.sh first (Section 10)." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: DIAMOND blastx vs CAZy (all kingdoms), R1 and R2 separately,
# then merge. staxids (col 13) carries the NCBI taxid of the best-hit
# sequence for kingdom assignment in 08_integrate.R.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_HITS}" && -s "${OUT_HITS}" ]]; then
    echo "[SKIP] CAZy readmap output already exists"
else
    if [[ -f "${OUT_R1}" && -s "${OUT_R1}" ]]; then
        echo "--- Step 1a: [SKIP] R1 output already exists ---"
    else
        echo "--- Step 1a: DIAMOND blastx vs CAZy (R1) ---"
        diamond blastx \
            --query "${R1_FASTQ}" \
            --db "${CAZY_DB}" \
            --out "${OUT_R1}" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                       qstart qend sstart send evalue bitscore staxids \
            --evalue ${EVALUE} \
            --max-target-seqs 1 \
            --sensitive \
            --threads ${THREADS}
        echo "Step 1a done: $(date)"
    fi

    if [[ -f "${OUT_R2}" && -s "${OUT_R2}" ]]; then
        echo "--- Step 1b: [SKIP] R2 output already exists ---"
    else
        echo "--- Step 1b: DIAMOND blastx vs CAZy (R2) ---"
        diamond blastx \
            --query "${R2_FASTQ}" \
            --db "${CAZY_DB}" \
            --out "${OUT_R2}" \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen \
                       qstart qend sstart send evalue bitscore staxids \
            --evalue ${EVALUE} \
            --max-target-seqs 1 \
            --sensitive \
            --threads ${THREADS}
        echo "Step 1b done: $(date)"
    fi

    # Ensure a trailing newline after R1 before appending R2 — DIAMOND does
    # not guarantee a terminal newline; bare cat R1 R2 can merge the last R1
    # line with the first R2 line into one malformed record.
    { cat "${OUT_R1}"; echo; cat "${OUT_R2}"; } > "${OUT_HITS}"
    rm -f "${OUT_R1}" "${OUT_R2}"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# ─────────────────────────────────────────────────────────────────────────────

N_HITS=$(wc -l < "${OUT_HITS}")
N_UNIQUE_CAZY=$(cut -f2 "${OUT_HITS}" | sort -u | wc -l)
N_UNIQUE_TAXA=$(cut -f13 "${OUT_HITS}" | sort -u | wc -l)

echo "  Total read hits (R1+R2)        : ${N_HITS}"
echo "  Unique CAZy sequences hit       : ${N_UNIQUE_CAZY}"
echo "  Unique taxids in hits           : ${N_UNIQUE_TAXA}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OUT_HITS}"
echo "============================================================"
