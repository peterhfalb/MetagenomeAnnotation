#!/bin/bash
# =============================================================================
# 04_dbcan.sh
# Slurm array job: CAZyme annotation with dbCAN3.
#
# dbCAN3 identifies carbohydrate-active enzymes (CAZymes) by running three
# evidence modes:
#   HMMER     — searched against dbCAN HMM profiles (CAZy family-level)
#   DIAMOND   — searched against characterized CAZy proteins
#   dbCAN_sub — searched against substrate-annotated CAZy subfamily HMMs;
#               gives a per-protein predicted glycan substrate (e.g. cellulose,
#               xylan, chitin) independent of genomic context, so it works for
#               fungal genomes (unlike dbCAN3's other substrate method — CGC/PUL
#               gene-cluster homology — which models bacterial-style
#               Polysaccharide Utilization Loci and isn't used here).
# Genes called by ≥2 of the family-calling tools (HMMER, DIAMOND; field
# #ofTools ≥ 2) are considered reliable family assignments. dbCAN_sub's
# substrate call is reported separately and doesn't affect that count.
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
#   overview.txt   main results: gene_id, HMMER call, DIAMOND call, #ofTools,
#                  dbCAN_sub substrate call (column name/presence depends on
#                  the installed run_dbcan version — see NOTE below)
#   hmmer.out      raw HMMER output
#   diamond.out    raw DIAMOND output
#
# NOTE: the exact --methods value for the dbCAN_sub evidence stream
# ("dbCANsub" below) should be verified against the installed run_dbcan
# version with `run_dbcan CAZyme_annotation --help` before relying on this in
# production — dbCAN's CLI flag naming has shifted across releases.
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

# dbCAN v4 changed the overview filename across sub-versions; detect it.
# Searches up to two directory levels deep for any known candidate.
find_overview() {
    find "${OUT_DIR}" -maxdepth 2 \
        \( -name "overview.txt" -o -name "overview.tsv" -o -name "CAZyme_annotation.tsv" \) \
        -size +0c 2>/dev/null | head -1
}

OVERVIEW=""
OVERVIEW=$(find_overview) || true

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
if [[ ! -f "${DB_DIR}/dbcan/dbCAN_sub.hmm.h3i" ]]; then
    echo "ERROR: dbCAN_sub database not found in ${DB_DIR}/dbcan/" >&2
    echo "  Run 00_setup_databases.sh first (Section 4 now presses dbCAN_sub.hmm)." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: run_dbcan
# Input type 'protein' is appropriate for eukaryotic gene predictions from
# MetaEuk (as opposed to 'prok' which skips signal peptide prediction, or
# 'meta' which runs additional steps for mixed-domain metagenomes).
# --methods hmm,diamond,dbCANsub: family calling (hmm, diamond) plus the
#   substrate-prediction evidence stream (dbCANsub). SignalP excluded
#   (requires a separate license and is less critical for fungal CAZyme
#   calling). CGC/PUL substrate prediction (gene-cluster homology to
#   bacterial-style PULs) is intentionally NOT run — not applicable to
#   fungal genomes, and would need the TF.hmm/STP.hmm/dbCAN-PUL.xlsx files
#   that 00_setup_databases.sh deliberately stubs out.
# ─────────────────────────────────────────────────────────────────────────────

if [[ -n "${OVERVIEW}" ]]; then
    echo "[SKIP] dbCAN3 output already exists: $(basename "${OVERVIEW}")"
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
        --methods hmm,diamond,dbCANsub \
        --threads ${THREADS}
    echo "Step 1 done: $(date)"

    OVERVIEW=$(find_overview) || true
    if [[ -z "${OVERVIEW}" ]]; then
        echo "ERROR: run_dbcan completed but no overview file found in ${OUT_DIR}" >&2
        echo "  Files present:" >&2
        ls -la "${OUT_DIR}" >&2
        exit 1
    fi
    echo "  Overview file: $(basename "${OVERVIEW}")"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# High-confidence CAZymes = supported by ≥2 tools (#ofTools column ≥ 2).
# Substrate count is purely informational here — exact column name/position
# depends on the installed run_dbcan version (see NOTE at top of file);
# 08_integrate.R locates it dynamically by name pattern rather than position.
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: CAZyme summary ---"

N_TOTAL=$(grep -c '>' "${PROTEINS_FAA}" 2>/dev/null || echo 0)
# overview.txt has a header line; skip it with NR>1
N_CALLED=$(awk 'NR>1 {c++} END {print c+0}' "${OVERVIEW}")
N_HICOF=$(awk 'NR>1 && $NF >= 2 {c++} END {print c+0}' "${OVERVIEW}")
N_SUBSTRATE=$(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if (tolower($i) ~ /substrate/) col=i}
                          NR>1 && col && $col != "-" && $col != "" {c++} END {print c+0}' "${OVERVIEW}")

echo "  Proteins searched           : ${N_TOTAL}"
echo "  Total CAZyme calls          : ${N_CALLED}"
echo "  High-confidence (≥2 tools)  : ${N_HICOF}"
echo "  Substrate calls (dbCAN_sub) : ${N_SUBSTRATE}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key output : ${OVERVIEW}"
echo "============================================================"
