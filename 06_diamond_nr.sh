#!/bin/bash
# =============================================================================
# 06_diamond_nr.sh
# Slurm array job: protein-level taxonomic classification via DIAMOND vs NR.
#
# Searches MetaEuk-predicted proteins against the NCBI NR protein database.
# Returns the top 10 hits per protein with taxonomic metadata. The taxonomy
# is embedded in the DIAMOND database (built with --taxonmap in step 00).
#
# Output columns (DIAMOND outfmt 6):
#   qseqid    query protein ID (MetaEuk gene ID)
#   sseqid    NR accession of the subject
#   pident    percent identity
#   length    alignment length
#   evalue    E-value
#   bitscore
#   staxids   NCBI taxon ID(s) of the subject (semicolon-separated if multiple)
#   sscinames scientific name(s) of the subject
#   sskingdoms kingdom classification (Eukaryota / Bacteria / Archaea / Viruses)
#   stitle    subject sequence title (species, gene function from NR header)
#
# Downstream use (in R):
#   - Filter hits where sskingdoms == "Eukaryota" to confirm eukaryotic origin
#   - Use staxids with the NCBI taxonomy tree to assign LCA at genus/family/order
#   - stitle provides quick human-readable identification of the top hit
#
# NOTE: This is the slowest annotation step (~8-16 hr per sample).
# It runs at --sensitive mode for better recovery of divergent fungal proteins.
# --evalue 1e-3 is intentionally lenient; stringent filtering happens in R
# during LCA calculation where only hits within a bitscore fraction of the
# best hit are included.
#
# Do not run directly. Submit via submit_diamond_nr.sh.
# =============================================================================

#SBATCH --job-name=metaG_diamond_nr
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
module load diamond/2.0.15-gcc-8.2.0-gkldzx7

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=32
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

# Top N hits per query for LCA calculation. 10 is standard — enough to
# resolve the LCA without inflating runtime with distant/uninformative hits.
MAX_TARGETS=10

# Lenient E-value for initial retrieval; stringent filtering happens in R.
EVALUE=1e-3

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

PROTEINS_FAA="${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas"
OUT_TSV="${ANNOTATION_DIR}/diamond_nr/${SAMPLE}_diamond_nr.tsv"
TMP_DIR="/scratch.global/${USER}/metaG_diamond_tmp/${SAMPLE}"

echo "============================================================"
echo "Sample     : ${SAMPLE}"
echo "Job ID     : ${SLURM_JOB_ID}  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start      : $(date)"
echo "============================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "${PROTEINS_FAA}" ]]; then
    echo "ERROR: MetaEuk proteins not found: ${PROTEINS_FAA}" >&2
    exit 1
fi
if [[ ! -s "${PROTEINS_FAA}" ]]; then
    echo "WARNING: Protein file is empty for ${SAMPLE} — skipping DIAMOND NR search."
    touch "${OUT_TSV}"
    exit 0
fi
if [[ ! -f "${DB_DIR}/diamond/nr.dmnd" ]]; then
    echo "ERROR: DIAMOND NR database not found: ${DB_DIR}/diamond/nr.dmnd" >&2
    echo "  Run 00_setup_databases.sh first." >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_TSV}")" "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: DIAMOND blastp vs NCBI NR
# --sensitive: required for soil fungi, many of which are novel or poorly
#   represented in NR. Without sensitivity mode, divergent fungal proteins
#   may receive no hits, leaving them taxonomically unclassified.
# --max-target-seqs 10: top 10 hits per protein for LCA calculation.
# staxids/sscinames/sskingdoms: taxonomy from the embedded NR taxon map.
# --tmpdir: DIAMOND writes temp files here; keep off of the network FS.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_TSV}" && -s "${OUT_TSV}" ]]; then
    echo "[SKIP] DIAMOND NR output already exists"
else
    echo "--- Step 1: DIAMOND blastp vs NR (~8-16 hr) ---"
    diamond blastp \
        --query "${PROTEINS_FAA}" \
        --db "${DB_DIR}/diamond/nr.dmnd" \
        --out "${OUT_TSV}" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore staxids sscinames sskingdoms stitle \
        --evalue ${EVALUE} \
        --max-target-seqs ${MAX_TARGETS} \
        --sensitive \
        --threads ${THREADS} \
        --tmpdir "${TMP_DIR}"
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# ─────────────────────────────────────────────────────────────────────────────

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
N_QUERIES_HIT=$(cut -f1 "${OUT_TSV}" | sort -u | wc -l)
N_EUK_HITS=$(awk '$9 == "Eukaryota" {c++} END {print c+0}' "${OUT_TSV}")

echo "  Proteins searched           : ${N_TOTAL}"
echo "  Proteins with ≥1 NR hit     : ${N_QUERIES_HIT}"
echo "  Hits with Eukaryota kingdom : ${N_EUK_HITS}"

rm -rf "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OUT_TSV}"
echo "============================================================"
