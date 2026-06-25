#!/bin/bash
# =============================================================================
# 04_dbcan.sh
# Slurm array job: CAZyme annotation with dbCAN (standalone `dbcan` package,
# confirmed installed version 5.2.9 — a ground-up rewrite of the older v3/v4
# tool this script originally targeted; output file layout is different and
# is what this script's logic is now written against).
#
# Runs three evidence modes via `--methods hmm,diamond,dbCANsub`:
#   hmm       — HMMER vs dbCAN HMM profiles (CAZy family-level)
#   diamond   — DIAMOND vs characterized CAZy proteins (CAZy family-level)
#   dbCANsub  — HMMER vs substrate-annotated CAZy subfamily HMMs; gives a
#               per-protein predicted glycan substrate (e.g. cellulose, xylan,
#               chitin) independent of genomic context, so it works for
#               fungal genomes (unlike dbCAN's other substrate method —
#               CGC/PUL gene-cluster homology — which models bacterial-style
#               Polysaccharide Utilization Loci and isn't used here).
#
# Unlike the older dbCAN versions, this version does NOT write a unified
# overview.txt — each method writes its own file, and there is no built-in
# "#ofTools" agreement count. We compute that ourselves in Step 2 below (and
# again, per-gene, in 08_integrate.R) by checking which gene IDs appear in
# both the hmm and diamond result files.
#
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
# Outputs (in ${ANNOTATION_DIR}/dbcan/${SAMPLE}/), all tab-separated with a
# header row, "Target Name"/"Gene ID" containing the full MetaEuk FASTA
# header (normalized to gene_id by 08_integrate.R):
#   dbCAN_hmm_results.tsv    HMM Name, HMM Length, Target Name, Target Length,
#                            i-Evalue, HMM From, HMM To, Target From, Target
#                            To, Coverage, HMM File Name
#   diamond.out              Gene ID, CAZy ID, % Identical, Length, ...
#                            (CAZy ID embeds family as "ACCESSION|FAMILY")
#   dbCANsub_hmm_results.tsv Subfam Name, Subfam Composition, Subfam EC,
#                            Substrate, HMM Length, Target Name, ... — the
#                            substrate-mapped result (via fam-substrate-
#                            mapping.tsv); use this file, not _raw.tsv.
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

# ── Conda (dbCAN is in metaG_annotation environment) ──────────────────────────
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

HMM_OUT="${OUT_DIR}/dbCAN_hmm_results.tsv"
DIAMOND_OUT="${OUT_DIR}/diamond.out"
SUB_OUT="${OUT_DIR}/dbCANsub_hmm_results.tsv"

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
    touch "${HMM_OUT}" "${DIAMOND_OUT}" "${SUB_OUT}"
    exit 0
fi
if [[ ! -f "${DB_DIR}/dbcan/dbCAN.hmm.h3i" ]]; then
    echo "ERROR: dbCAN database not found in ${DB_DIR}/dbcan/" >&2
    echo "  Run 00_setup_databases.sh first." >&2
    exit 1
fi
if [[ ! -f "${DB_DIR}/dbcan/dbCAN-sub.hmm.h3i" ]]; then
    echo "ERROR: dbCAN-sub database not found in ${DB_DIR}/dbcan/" >&2
    echo "  Run 00_setup_databases.sh first (Section 4 presses dbCAN-sub.hmm)." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: run_dbcan
# Input type 'protein' is appropriate for eukaryotic gene predictions from
# MetaEuk (as opposed to 'prok' which skips signal peptide prediction, or
# 'meta' which runs additional steps for mixed-domain metagenomes).
# CGC/PUL substrate prediction (gene-cluster homology to bacterial-style
# PULs) is intentionally NOT requested — not applicable to fungal genomes.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${HMM_OUT}" && -f "${DIAMOND_OUT}" && -f "${SUB_OUT}" ]]; then
    echo "[SKIP] dbCAN output already exists"
else
    echo "--- Step 1: run_dbcan ---"
    run_dbcan CAZyme_annotation \
        --input_raw_data "${PROTEINS_FAA}" \
        --mode protein \
        --db_dir "${DB_DIR}/dbcan" \
        --output_dir "${OUT_DIR}" \
        --methods hmm,diamond,dbCANsub \
        --threads ${THREADS}
    echo "Step 1 done: $(date)"

    if [[ ! -f "${HMM_OUT}" || ! -f "${DIAMOND_OUT}" || ! -f "${SUB_OUT}" ]]; then
        echo "ERROR: run_dbcan completed but expected output files are missing in ${OUT_DIR}" >&2
        echo "  Expected: $(basename "${HMM_OUT}"), $(basename "${DIAMOND_OUT}"), $(basename "${SUB_OUT}")" >&2
        echo "  Files present:" >&2
        ls -la "${OUT_DIR}" >&2
        exit 1
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# No unified "#ofTools" column in this dbCAN version — compute the HMMER ∩
# DIAMOND agreement ourselves by comparing gene IDs (column 3 of the hmm
# file, column 1 of the diamond file) between the two result files.
# Substrate count is purely informational here; 08_integrate.R does the real
# per-gene join.
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: CAZyme summary ---"

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)

HMM_IDS=$(awk -F'\t' 'NR>1 {print $3}' "${HMM_OUT}" | sort -u)
DIAMOND_IDS=$(awk -F'\t' 'NR>1 {print $1}' "${DIAMOND_OUT}" | sort -u)

N_HMM=$(printf '%s\n' "${HMM_IDS}" | grep -c . || true)
N_DIAMOND=$(printf '%s\n' "${DIAMOND_IDS}" | grep -c . || true)
N_HICOF=$(comm -12 <(printf '%s\n' "${HMM_IDS}") <(printf '%s\n' "${DIAMOND_IDS}") | grep -c . || true)
N_SUBSTRATE=$(awk -F'\t' 'NR>1 && $4 != "" && $4 != "-" {c++} END {print c+0}' "${SUB_OUT}")

echo "  Proteins searched           : ${N_TOTAL}"
echo "  HMMER family calls          : ${N_HMM}"
echo "  DIAMOND family calls        : ${N_DIAMOND}"
echo "  High-confidence (HMM∩DIAMOND): ${N_HICOF}"
echo "  Substrate calls (dbCAN_sub)  : ${N_SUBSTRATE}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key outputs:"
echo "  ${HMM_OUT}"
echo "  ${DIAMOND_OUT}"
echo "  ${SUB_OUT}"
echo "============================================================"
