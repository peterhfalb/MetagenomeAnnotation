#!/bin/bash
# =============================================================================
# 06c_mycocosm_taxonomy.sh
# Slurm array job: protein-level taxonomic classification via MMseqs2 against
# a custom Mycocosm (fungal) + Phytozome (plant) reference database.
#
# Runs ADDITIONALLY alongside 06_mmseqs_taxonomy.sh (UniRef90) — this step
# does not replace it. UniRef90 gives broad kingdom-level taxonomy/contamination
# checks; this step gives finer fungal subphylum resolution (Agaricomycotina,
# Pezizomycotina) and flags plant-derived sequences, replicating the taxonomy
# approach used to identify EM-forming fungi in the source methodology paper.
#
# NO FILTERING is applied anywhere based on this taxonomy — it only adds
# resolution. All genes are retained regardless of their classification here.
#
# Same easy-taxonomy mechanics as 06_mmseqs_taxonomy.sh (2bLCA, --tax-lineage),
# just pointed at the custom DB built in 00_setup_databases.sh Section 8.
#
# Inputs:
#   ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas   (MetaEuk proteins)
#   ${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome/db     (built in step 00)
#
# Outputs:
#   ${ANNOTATION_DIR}/mycocosm_taxonomy/${SAMPLE}_lca.tsv
#     Columns: gene_id, lca_taxid, lca_rank, lca_name, lineage
#     Identical schema to 06_mmseqs_taxonomy.sh's output, so the same R
#     reader shape can be reused (with a myco_ column prefix downstream).
#
# Do not run directly. Submit via submit_mycocosm_taxonomy.sh.
# =============================================================================

#SBATCH --job-name=metaG_myco_tax
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

# NOTE on resources: this DB is expected to be far smaller than UniRef90
# (~30 GB), since it's a curated panel of JGI genomes rather than all of
# UniRef90. The 128 GB / --split-memory-limit settings below are inherited
# conservatively from 06_mmseqs_taxonomy.sh — once the real DB size is known,
# both can likely be reduced.

set -euo pipefail

# ── Environment ───────────────────────────────────────────────────────────────
set +u
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_mmseqs
set -u

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=32
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

SENSITIVITY=6
LCA_MODE=2

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

PROTEINS_FAA="${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas"
OUT_DIR="${ANNOTATION_DIR}/mycocosm_taxonomy"
OUT_PREFIX="${OUT_DIR}/${SAMPLE}"
OUT_LCA="${OUT_PREFIX}_lca.tsv"
TMP_DIR="/scratch.global/${USER}/metaG_myco_tmp/${SAMPLE}"

MYCO_DB="${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome/db"

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
if [[ ! -f "${MYCO_DB}" ]]; then
    echo "ERROR: Mycocosm/Phytozome MMseqs2 database not found: ${MYCO_DB}" >&2
    echo "  Requires a manual download + taxid mapping step — see SECTION 8 in" >&2
    echo "  00_setup_databases.sh for instructions." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: mmseqs easy-taxonomy vs the custom Mycocosm/Phytozome database.
# Same flags as 06_mmseqs_taxonomy.sh: --lca-mode 2 (2bLCA), --tax-lineage 1
# (full semicolon-separated lineage), -s 6 (sensitive). lca_taxid 0 = no hit.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_LCA}" && -s "${OUT_LCA}" ]]; then
    echo "[SKIP] Mycocosm/Phytozome taxonomy output already exists"
else
    echo "--- Step 1: mmseqs easy-taxonomy (Mycocosm/Phytozome) ---"
    mmseqs easy-taxonomy \
        "${PROTEINS_FAA}" \
        "${MYCO_DB}" \
        "${OUT_PREFIX}" \
        "${TMP_DIR}" \
        --threads ${THREADS} \
        --lca-mode ${LCA_MODE} \
        --tax-lineage 1 \
        -s ${SENSITIVITY} \
        --split-memory-limit 100G
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary — purely informational, no filtering is applied.
# ─────────────────────────────────────────────────────────────────────────────

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
N_CLASSIFIED=$(awk -F'\t' '$2 != "0" && $2 != "" {c++} END {print c+0}' "${OUT_LCA}")
N_PLANT=$(awk -F'\t' '$5 ~ /Viridiplantae/ {c++} END {print c+0}' "${OUT_LCA}")
N_FUNGAL_SUBPHYLUM=$(awk -F'\t' '$5 ~ /Agaricomycotina|Pezizomycotina/ {c++} END {print c+0}' "${OUT_LCA}")

echo "  Proteins searched              : ${N_TOTAL}"
echo "  Classified (taxid != 0)        : ${N_CLASSIFIED}"
echo "  Plant (Viridiplantae) in lineage: ${N_PLANT}"
echo "  Agaricomycotina/Pezizomycotina  : ${N_FUNGAL_SUBPHYLUM}"

rm -rf "${TMP_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OUT_LCA}"
echo "============================================================"
