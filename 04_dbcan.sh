#!/bin/bash
# =============================================================================
# 04_dbcan.sh
# Slurm array job: CAZyme annotation with dbCAN3.
#
# dbCAN3 identifies carbohydrate-active enzymes (CAZymes) by running two
# evidence modes and reporting genes supported by both as high-confidence:
#   HMMER  — searched against dbCAN HMM profiles (CAZy family-level)
#   DIAMOND — searched against characterized CAZy proteins
# Genes called by ≥2 tools (field #ofTools ≥ 2) are considered reliable.
# CAZy families relevant to carbon acquisition in fungi:
#   GH (glycoside hydrolases)   — cellulose, hemicellulose, starch degradation
#   PL (polysaccharide lyases)  — pectin degradation
#   CE (carbohydrate esterases) — deacetylation of plant polymers
#   AA (auxiliary activities)   — oxidative degradation (lignin, cellulose)
#   CBM (binding modules)       — substrate targeting, not catalytic
#
# Inputs:
#   ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas   (MetaEuk proteins)
#
# Outputs (in ${ANNOTATION_DIR}/dbcan/${SAMPLE}/):
#   overview.txt   main results: gene_id, HMMER call, DIAMOND call, #ofTools
#   hmmer.out      raw HMMER output
#   diamond.out    raw DIAMOND output
#
# Do not run directly. Submit via submit_dbcan.sh.
# =============================================================================

#SBATCH --job-name=metaG_dbcan
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Conda (dbCAN3 is in metaG_annotation environment) ────────────────────────
set +u
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_annotation
set -u

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=16
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

PROTEINS_FAA="${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.fas"
OUT_DIR="${ANNOTATION_DIR}/dbcan/${SAMPLE}"
OVERVIEW="${OUT_DIR}/overview.txt"

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
    echo "WARNING: Protein file is empty for ${SAMPLE} — skipping dbCAN."
    mkdir -p "${OUT_DIR}"
    touch "${OUT_DIR}/overview.txt"
    exit 0
fi
if [[ ! -f "${DB_DIR}/dbcan/dbCAN.hmm.h3i" ]]; then
    echo "ERROR: dbCAN3 database not found in ${DB_DIR}/dbcan/" >&2
    echo "  Run 00_setup_databases.sh first." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: run_dbcan
# Input type 'protein' is appropriate for eukaryotic gene predictions from
# MetaEuk (as opposed to 'prok' which skips signal peptide prediction, or
# 'meta' which runs additional steps for mixed-domain metagenomes).
# --tools hmmer diamond: run both evidence modes; SignalP excluded (requires
#   a separate license and is less critical for fungal CAZyme calling).
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OVERVIEW}" && -s "${OVERVIEW}" ]]; then
    echo "[SKIP] dbCAN3 output already exists"
else
    echo "--- Step 1: run_dbcan ---"
    # dbCAN v4 replaced the old positional-argument syntax with subcommands.
    # Old (v3): run_dbcan <input> protein --db_dir ...
    # New (v4): run_dbcan CAZyme_annotation --input_raw_data <input> --mode protein ...
    run_dbcan CAZyme_annotation \
        --input_raw_data "${PROTEINS_FAA}" \
        --mode protein \
        --db_dir "${DB_DIR}/dbcan" \
        --output_dir "${OUT_DIR}" \
        --methods hmm,diamond \
        --threads ${THREADS}
    echo "Step 1 done: $(date)"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# High-confidence CAZymes = supported by ≥2 tools (#ofTools column ≥ 2).
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: CAZyme summary ---"

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
# overview.txt has a header line; skip it with NR>1
N_CALLED=$(awk 'NR>1 {c++} END {print c+0}' "${OVERVIEW}")
N_HICOF=$(awk 'NR>1 && $NF >= 2 {c++} END {print c+0}' "${OVERVIEW}")

echo "  Proteins searched           : ${N_TOTAL}"
echo "  Total CAZyme calls          : ${N_CALLED}"
echo "  High-confidence (≥2 tools)  : ${N_HICOF}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OVERVIEW}"
echo "============================================================"
