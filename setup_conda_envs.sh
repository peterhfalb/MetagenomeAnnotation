#!/bin/bash
# =============================================================================
# setup_conda_envs.sh
# Creates conda environments for the metagenome annotation pipeline.
#
# Run ONCE interactively on the login node BEFORE submitting any jobs.
# Do NOT submit this to Slurm — conda environment creation needs login-node
# network access and should not be run as a batch job.
#
# Usage:
#   bash setup_conda_envs.sh
#
# Creates three environments:
#   metaG_tiara      — contig-level taxonomic classification (PyTorch dependency
#                      isolated to avoid conflicts with bioconda tools)
#   metaG_annotation — all other annotation tools: KOfamScan, dbCAN3,
#                      featureCounts (subread), MMseqs2 (CPU, for DB builds), HMMER
#   metaG_mmseqs     — GPU-enabled MMseqs2 for taxonomy (CUDA-compiled; isolated
#                      because CUDA runtime dependencies can conflict with other
#                      bioconda packages). At runtime, load cuda/12.1.1 alongside.
#
# Tools loaded as MSI modules (NOT installed via conda):
#   MetaEuk  — metaeuk/6-a5d39d9-gcc-8.2.0-ji6jath
#   DIAMOND  — diamond/2.0.15-gcc-8.2.0-gkldzx7
#   samtools — samtools/1.21
#   CUDA     — cuda/12.1.1  (loaded at runtime for metaG_mmseqs jobs)
# =============================================================================

set -euo pipefail

source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate base

# ── Tiara environment ─────────────────────────────────────────────────────────
# Tiara uses PyTorch for deep-learning contig classification. Isolated to
# prevent torch solver conflicts with bioconda packages.

if conda env list | grep -qx "metaG_tiara .*"; then
    echo "[SKIP] metaG_tiara already exists"
else
    echo "--- Creating metaG_tiara environment ---"
    mamba create -y -n metaG_tiara \
        -c conda-forge -c bioconda \
        tiara
    echo "metaG_tiara done"
fi

# ── Main annotation environment ───────────────────────────────────────────────
# kofamscan  : KEGG KO assignment via HMM profiles
# dbcan      : CAZyme annotation (run_dbcan, includes HMMER and DIAMOND modes)
# subread    : provides featureCounts for read-to-gene quantification
# mmseqs2    : required to build the MetaEuk fungal reference database
# hmmer      : HMMER3 for KOfamScan and dbCAN HMM searches

if conda env list | grep -qx "metaG_annotation .*"; then
    echo "[SKIP] metaG_annotation already exists"
else
    echo "--- Creating metaG_annotation environment ---"
    mamba create -y -n metaG_annotation \
        -c conda-forge -c bioconda \
        kofamscan \
        dbcan \
        subread \
        mmseqs2 \
        hmmer
    echo "metaG_annotation done"
fi

# ── MMseqs2 GPU environment ───────────────────────────────────────────────────
# GPU-accelerated MMseqs2 for the taxonomy step. The conda package includes
# CUDA-compiled binaries; load cuda/12.1.1 at runtime so the dynamic linker
# finds the CUDA libraries. Isolated from metaG_annotation to avoid conflicts.

if conda env list | grep -qx "metaG_mmseqs .*"; then
    echo "[SKIP] metaG_mmseqs already exists"
else
    echo "--- Creating metaG_mmseqs environment ---"
    mamba create -y -n metaG_mmseqs \
        -c conda-forge -c bioconda \
        mmseqs2
    echo "metaG_mmseqs done"
fi

echo ""
echo "============================================================"
echo "Environments ready."
echo "  Activate Tiara  : conda activate metaG_tiara"
echo "  Activate main   : conda activate metaG_annotation"
echo "  Activate MMseqs2: conda activate metaG_mmseqs"
echo "  (load cuda/12.1.1 alongside metaG_mmseqs for GPU jobs)"
echo "============================================================"
echo ""
echo "Next: submit database setup job with:"
echo "  bash submit_setup_databases.sh"
