#!/bin/bash
# =============================================================================
# 06_mmseqs_taxonomy.sh
# Slurm array job: protein-level taxonomic classification via MMseqs2.
#
# Uses mmseqs easy-taxonomy to search MetaEuk-predicted proteins against
# UniRef90 and assign taxonomy via 2bLCA (two-way best-hit LCA). The LCA is
# computed internally by MMseqs2 — no per-hit post-processing in R required.
# GPU acceleration is used if cuda/12.1.1 and a GPU node are available.
#
# Inputs:
#   ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas   (MetaEuk proteins)
#   ${DB_DIR}/mmseqs_taxonomy/uniref90                  (built in step 00)
#
# Outputs:
#   ${ANNOTATION_DIR}/mmseqs_taxonomy/${SAMPLE}_lca.tsv
#     Columns: gene_id, lca_taxid, lca_rank, lca_name, lineage
#     One row per predicted protein. lca_taxid = 0 means unclassified.
#     The lineage column (--tax-lineage 1) gives the full path as
#     semicolon-separated taxon names (root to LCA).
#
# Do not run directly. Submit via submit_mmseqs_taxonomy.sh.
# =============================================================================

#SBATCH --job-name=metaG_mmseqs_tax
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

# NOTE on GPU: MMseqs2 v18 uses a gpuserver client-server architecture rather
# than a --gpu flag on individual commands. Running a persistent GPU server per
# SLURM array task is non-trivial. CPU MMseqs2 is already 10-20x faster than
# DIAMOND for taxonomy — this step should complete in 30-90 min/sample on 32
# CPUs. GPU acceleration can be revisited later if throughput is a bottleneck.

set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────────
set +u
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_mmseqs
set -u

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=16
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

# Sensitivity: 1 (fast) to 7.5 (max). 6 is a good metagenomics default.
# With GPU this runs fast even at higher sensitivity — increase to 7.5 if
# you want maximum recovery of divergent fungal proteins.
SENSITIVITY=6

# LCA mode: 2 = 2bLCA (top hit + LCA). Recommended for metagenomics —
# more specific than pure LCA while still robust to database contamination.
LCA_MODE=2

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

PROTEINS_FAA="${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas"
OUT_DIR="${ANNOTATION_DIR}/mmseqs_taxonomy"
OUT_PREFIX="${OUT_DIR}/${SAMPLE}"
OUT_LCA="${OUT_PREFIX}_lca.tsv"
TMP_DIR="/scratch.global/${USER}/metaG_mmseqs_tmp/${SAMPLE}"

echo "============================================================"
echo "Sample     : ${SAMPLE}"
echo "Job ID     : ${SLURM_JOB_ID}  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start      : $(date)"
echo "MMseqs2    : $(mmseqs version)"
echo "============================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "${PROTEINS_FAA}" ]]; then
    echo "ERROR: MetaEuk proteins not found: ${PROTEINS_FAA}" >&2
    exit 1
fi
if [[ ! -s "${PROTEINS_FAA}" ]]; then
    echo "WARNING: Protein file is empty for ${SAMPLE} — skipping taxonomy."
    touch "${OUT_LCA}"
    exit 0
fi
if [[ ! -f "${DB_DIR}/mmseqs_taxonomy/uniref90" ]]; then
    echo "ERROR: MMseqs2 UniRef90 database not found: ${DB_DIR}/mmseqs_taxonomy/uniref90" >&2
    echo "  Run 00_setup_databases.sh first." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: mmseqs easy-taxonomy
# easy-taxonomy runs the full pipeline internally:
#   createdb → search → lca → convertalis
# --lca-mode 2    : 2bLCA — combines top-hit and LCA for a balanced assignment
# --tax-lineage 1 : adds full semicolon-separated lineage names to output
# -s              : sensitivity (6 = sensitive, good for metagenomes)
#
# Output ${SAMPLE}_lca.tsv columns:
#   query_id  lca_taxid  lca_rank  lca_name  lineage
# lca_taxid 0 = unclassified (no hits above threshold).
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_LCA}" && -s "${OUT_LCA}" ]]; then
    echo "[SKIP] MMseqs2 taxonomy output already exists"
else
    echo "--- Step 1: mmseqs easy-taxonomy ---"
    mmseqs easy-taxonomy \
        "${PROTEINS_FAA}" \
        "${DB_DIR}/mmseqs_taxonomy/uniref90" \
        "${OUT_PREFIX}" \
        "${TMP_DIR}" \
        --threads ${THREADS} \
        --lca-mode ${LCA_MODE} \
        --tax-lineage 1 \
        -s ${SENSITIVITY}
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# ─────────────────────────────────────────────────────────────────────────────

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
N_CLASSIFIED=$(awk -F'\t' '$2 != "0" && $2 != "" {c++} END {print c+0}' "${OUT_LCA}")
N_EUK=$(awk -F'\t' '$5 ~ /Eukaryota/ {c++} END {print c+0}' "${OUT_LCA}")
N_FUNGI=$(awk -F'\t' '$5 ~ /Fungi/ {c++} END {print c+0}' "${OUT_LCA}")

echo "  Proteins searched         : ${N_TOTAL}"
echo "  Classified (taxid != 0)   : ${N_CLASSIFIED}"
echo "  Eukaryota in lineage      : ${N_EUK}"
echo "  Fungi in lineage          : ${N_FUNGI}"

rm -rf "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OUT_LCA}"
echo "============================================================"
