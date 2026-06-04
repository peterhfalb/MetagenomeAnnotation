#!/bin/bash
# =============================================================================
# submit_setup_databases.sh
# Wrapper: validates prerequisites and submits 00_setup_databases.sh.
# Run once from the login node after setup_conda_envs.sh has completed.
#
# Usage:
#   bash submit_setup_databases.sh
#
# Prerequisites:
#   1. bash setup_conda_envs.sh         (creates conda environments)
#   2. PHI-base FASTA manually placed at:
#      /projects/standard/kennedyp/shared/databases/metaG_annotation/phibase/phi-base_current.fas
#      (see Section 5 in 00_setup_databases.sh for download instructions)
#   3. BLAST_MODULE in 00_setup_databases.sh updated with the correct module name
#      (find it with: module spider blast)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"
LOG_DIR="${SCRIPT_DIR}/logs/setup"

mkdir -p "${LOG_DIR}"

# ── Prerequisite checks ───────────────────────────────────────────────────────

source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate base

if ! conda env list | grep -qx "metaG_annotation .*"; then
    echo "ERROR: conda environment 'metaG_annotation' not found."
    echo "  Run first: bash setup_conda_envs.sh"
    exit 1
fi

if ! conda env list | grep -qx "metaG_tiara .*"; then
    echo "ERROR: conda environment 'metaG_tiara' not found."
    echo "  Run first: bash setup_conda_envs.sh"
    exit 1
fi

if [[ -f "${DB_DIR}/databases_complete.flag" ]]; then
    echo "WARNING: databases_complete.flag already exists — databases may already be built."
    echo "  Remove flag and re-run to force rebuild: rm ${DB_DIR}/databases_complete.flag"
    read -r -p "  Submit anyway? [y/N] " reply
    [[ "${reply}" =~ ^[Yy]$ ]] || exit 0
fi

# Check for PHI-base FASTA (warn but don't block)
PHIBASE_FASTA="${DB_DIR}/phibase/phi-base_current.fas"
if [[ ! -f "${PHIBASE_FASTA}" ]]; then
    echo ""
    echo "WARNING: PHI-base FASTA not found at:"
    echo "  ${PHIBASE_FASTA}"
    echo "  The job will still run but will skip the PHI-base DIAMOND build."
    echo "  See Section 5 in 00_setup_databases.sh for download instructions."
    echo "  You can re-run this submission after placing the FASTA to add PHI-base."
    echo ""
fi

# ── Submit ────────────────────────────────────────────────────────────────────

JOB_ID=$(sbatch \
    --output="${LOG_DIR}/setup_databases_%j.out" \
    --error="${LOG_DIR}/setup_databases_%j.err" \
    "${SCRIPT_DIR}/00_setup_databases.sh" \
    | awk '{print $NF}')

echo "Submitted database setup job: ${JOB_ID}"
echo "Logs: ${LOG_DIR}/setup_databases_${JOB_ID}.out"
echo ""
echo "Monitor with:"
echo "  squeue -j ${JOB_ID}"
echo "  tail -f ${LOG_DIR}/setup_databases_${JOB_ID}.out"
