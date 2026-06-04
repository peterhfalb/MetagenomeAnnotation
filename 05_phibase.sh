#!/bin/bash
# =============================================================================
# 05_phibase.sh
# Slurm array job: pathogenicity gene annotation via DIAMOND vs PHI-base.
#
# PHI-base (Pathogen-Host Interaction database) contains experimentally
# verified pathogenicity, virulence, and effector genes from fungal and
# oomycete pathogens. A DIAMOND search identifies predicted proteins with
# homology to known pathogenicity-related proteins.
#
# DIAMOND output columns (format 6):
#   qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
#
# The stitle field contains the PHI-base accession and phenotype annotation
# in the format: PHI:<accession>|gene_name|pathogen_species|host|phenotype|...
# This encodes whether the gene is required for pathogenicity, reduced virulence,
# loss of pathogenicity, etc. Parse stitle in R to extract phenotype class.
#
# Inputs:
#   ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas   (MetaEuk proteins)
#   ${DB_DIR}/phibase/phi-base.dmnd
#
# Outputs:
#   ${ANNOTATION_DIR}/phibase/${SAMPLE}_phibase.tsv   DIAMOND tabular hits
#
# Do not run directly. Submit via submit_phibase.sh.
# =============================================================================

#SBATCH --job-name=metaG_phibase
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Modules ───────────────────────────────────────────────────────────────────
module load diamond/2.0.15-gcc-8.2.0-gkldzx7

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=8
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

# E-value cutoff for PHI-base hits. 1e-5 is standard for homology-based
# pathogenicity annotation — stringent enough to avoid spurious hits while
# still capturing divergent effectors.
EVALUE=1e-5

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

PROTEINS_FAA="${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas"
OUT_TSV="${ANNOTATION_DIR}/phibase/${SAMPLE}_phibase.tsv"

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
    echo "WARNING: Protein file is empty for ${SAMPLE} — skipping PHI-base search."
    touch "${OUT_TSV}"
    exit 0
fi
if [[ ! -f "${DB_DIR}/phibase/phi-base.dmnd" ]]; then
    echo "ERROR: PHI-base DIAMOND database not found: ${DB_DIR}/phibase/phi-base.dmnd" >&2
    echo "  Complete the manual PHI-base download step (see 00_setup_databases.sh Section 5)." >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_TSV}")"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: DIAMOND blastp vs PHI-base
# --max-target-seqs 1: best PHI-base hit per query protein only.
#   Rationale: we want the most similar known pathogenicity gene for each
#   predicted protein; multiple hits add noise without improving classification.
# --sensitive: PHI-base is relatively small; sensitivity mode is fast here
#   and recovers more divergent pathogenicity gene homologs.
# stitle is included to capture the PHI-base phenotype annotation string.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_TSV}" && -s "${OUT_TSV}" ]]; then
    echo "[SKIP] PHI-base output already exists"
else
    echo "--- Step 1: DIAMOND blastp vs PHI-base ---"
    diamond blastp \
        --query "${PROTEINS_FAA}" \
        --db "${DB_DIR}/phibase/phi-base.dmnd" \
        --out "${OUT_TSV}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
        --evalue ${EVALUE} \
        --max-target-seqs 1 \
        --sensitive \
        --threads ${THREADS}
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# ─────────────────────────────────────────────────────────────────────────────

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
N_HITS=$(wc -l < "${OUT_TSV}")

echo "  Proteins searched  : ${N_TOTAL}"
echo "  PHI-base hits      : ${N_HITS} (E-value ≤ ${EVALUE})"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OUT_TSV}"
echo "============================================================"
