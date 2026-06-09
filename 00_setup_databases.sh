#!/bin/bash
# =============================================================================
# 00_setup_databases.sh
# Slurm job: download and index all annotation databases.
#
# Idempotent — re-running skips any step whose output already exists.
# Submit via submit_setup_databases.sh, not directly.
#
# Databases built here (all stored permanently in DB_DIR):
#   taxonomy/     NCBI taxdump (names.dmp, nodes.dmp) + prot.accession2taxid
#   diamond/      DIAMOND-format NR (~80 GB), built from MSI's BLAST NR copy
#   kofam/        KOfam HMM profiles + ko_list (KEGG KO assignment)
#   dbcan/        dbCAN3 CAZyme databases (HMM + DIAMOND modes)
#   phibase/      PHI-base DIAMOND database (pathogenicity genes)
#   metaeuk/      Fungi RefSeq proteins + MMseqs2 database (MetaEuk gene calling)
#
# NOTE — PHI-base requires a manual download step before this job runs.
#        See SECTION 5 below for instructions.
#
# NOTE — DIAMOND NR build requires the BLAST+ module for blastdbcmd.
#        Find the correct module name: module spider blast
#        Then update BLAST_MODULE below.
# =============================================================================

#SBATCH --job-name=metaG_setup_dbs
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256gb
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=32

# Permanent database storage (shared lab folder)
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

# Scratch for large intermediate files (auto-cleared every 30 days — fine here)
SCRATCH_DIR="/scratch.global/falb0011"

# ── Environment ───────────────────────────────────────────────────────────────
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_annotation

# DIAMOND is still needed for the PHI-base and fungi DIAMOND database builds.
# BLAST+ module no longer required (replaced by mmseqs databases for taxonomy).
module load diamond/2.0.15-gcc-8.2.0-gkldzx7

# ── Create directories ────────────────────────────────────────────────────────
mkdir -p "${DB_DIR}"/{taxonomy,mmseqs_taxonomy,diamond,kofam,dbcan,phibase,metaeuk}
mkdir -p "${SCRATCH_DIR}"/{fungi_refseq,mmseqs_dl_tmp}

echo "============================================================"
echo "Database setup started : $(date)"
echo "DB_DIR                 : ${DB_DIR}"
echo "SCRATCH_DIR            : ${SCRATCH_DIR}"
echo "============================================================"

# =============================================================================
# SECTION 1: NCBI Taxonomy
# names.dmp and nodes.dmp are used by the R integration script (taxonomizr)
# to expand MMseqs2 LCA taxon IDs to standard ranks (phylum, class, order...).
# prot.accession2taxid is NOT needed — that was only required for DIAMOND NR.
# =============================================================================

if [[ ! -f "${DB_DIR}/taxonomy/nodes.dmp" ]]; then
    echo "--- [1] Downloading NCBI taxdump ---"
    wget -q -O "${SCRATCH_DIR}/taxdump.tar.gz" \
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf "${SCRATCH_DIR}/taxdump.tar.gz" -C "${DB_DIR}/taxonomy/"
    rm "${SCRATCH_DIR}/taxdump.tar.gz"
    echo "taxdump done: $(date)"
else
    echo "[SKIP] taxonomy/nodes.dmp already exists"
fi

# =============================================================================
# SECTION 2: MMseqs2 UniRef90 taxonomy database
# Used by mmseqs easy-taxonomy (step 06) for protein-level LCA classification.
# mmseqs databases handles download + indexing automatically — no BLAST module
# or blastdbcmd required. UniRef90 (~30 GB indexed) provides comprehensive
# coverage of fungi, bacteria, archaea, plants, and animals at 90% identity
# clustering. MMseqs2 computes LCA internally; the R script (step 08) only
# needs to expand the single LCA taxon ID per gene to standard ranks.
# Runtime: ~2-4 hr download + indexing.
# =============================================================================

MMSEQS_TAX_DB="${DB_DIR}/mmseqs_taxonomy/uniref90"

if [[ ! -f "${MMSEQS_TAX_DB}" ]]; then
    echo "--- [2] Downloading and indexing MMseqs2 UniRef90 (~30 GB, 2-4 hr) ---"
    mmseqs databases UniRef90 \
        "${MMSEQS_TAX_DB}" \
        "${SCRATCH_DIR}/mmseqs_dl_tmp" \
        --threads ${THREADS}
    echo "MMseqs2 UniRef90 done: $(date)"
else
    echo "[SKIP] MMseqs2 UniRef90 database already exists"
fi

# =============================================================================
# SECTION 3: KOfam — KEGG ortholog HMM profiles
# Used by KOfamScan to assign KEGG KO numbers to predicted proteins.
# KO numbers map to KEGG pathways (carbon/nitrogen acquisition, etc.).
# No KEGG license required — profiles and ko_list are freely distributed
# by GenomeNet.
# =============================================================================

if [[ ! -d "${DB_DIR}/kofam/profiles" ]]; then
    echo "--- [3a] Downloading KOfam HMM profiles (~5 GB) ---"
    wget -q -O "${SCRATCH_DIR}/profiles.tar.gz" \
        ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
    tar -xzf "${SCRATCH_DIR}/profiles.tar.gz" -C "${DB_DIR}/kofam/"
    rm "${SCRATCH_DIR}/profiles.tar.gz"
    echo "KOfam profiles done: $(date)"
else
    echo "[SKIP] KOfam profiles already exist"
fi

if [[ ! -f "${DB_DIR}/kofam/ko_list" ]]; then
    echo "--- [3b] Downloading KOfam ko_list ---"
    wget -q -O "${DB_DIR}/kofam/ko_list.gz" \
        ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    gunzip "${DB_DIR}/kofam/ko_list.gz"
    echo "KOfam ko_list done: $(date)"
else
    echo "[SKIP] KOfam ko_list already exists"
fi

# =============================================================================
# SECTION 4: dbCAN3 — CAZyme databases
# run_dbcan operates in three modes simultaneously:
#   HMMER  — dbCAN HMM profiles (most sensitive, primary result)
#   DIAMOND — searched against characterized CAZy proteins (used as evidence)
# Genes are called CAZymes when ≥2 modes agree (dbCAN3 default).
# run_dbcan download_db handles downloading and indexing all required files.
# =============================================================================

if [[ ! -f "${DB_DIR}/dbcan/dbCAN.hmm.h3i" ]]; then
    echo "--- [4] Downloading dbCAN3 databases ---"
    # dbCAN 4.x renamed the subcommand from 'download_db' to 'database'
    run_dbcan database --db-dir "${DB_DIR}/dbcan"
    echo "dbCAN3 done: $(date)"
else
    echo "[SKIP] dbCAN3 databases already exist"
fi

# =============================================================================
# SECTION 5: PHI-base — pathogen-host interaction database
# Provides curated fungal pathogenicity genes for DIAMOND search.
# Used to identify genes with known roles in infection/virulence.
#
# *** MANUAL STEP REQUIRED BEFORE THIS JOB RUNS ***
# PHI-base does not allow direct automated download. Do the following once:
#   1. Go to: https://www.phi-base.org/
#   2. Download the current release FASTA (phi-base_current.fas or similar)
#   3. Copy it to MSI:
#        scp phi-base_current.fas agate.msi.umn.edu:${DB_DIR}/phibase/phi-base_current.fas
# This script will then build the DIAMOND database automatically.
# =============================================================================

PHIBASE_FASTA="${DB_DIR}/phibase/phi-base_current.fas"
PHIBASE_DMND="${DB_DIR}/phibase/phi-base.dmnd"

if [[ ! -f "${PHIBASE_FASTA}" ]]; then
    echo ""
    echo "WARNING: PHI-base FASTA not found — skipping PHI-base DIAMOND build."
    echo "  Download the current FASTA from https://www.phi-base.org/ and place at:"
    echo "  ${PHIBASE_FASTA}"
    echo "  Then re-run this script; the DIAMOND build will run automatically."
    echo ""
elif [[ ! -f "${PHIBASE_DMND}" ]]; then
    echo "--- [5] Building PHI-base DIAMOND database ---"
    diamond makedb \
        --in "${PHIBASE_FASTA}" \
        --db "${PHIBASE_DMND}" \
        --threads ${THREADS}
    echo "PHI-base DIAMOND done: $(date)"
else
    echo "[SKIP] PHI-base DIAMOND database already exists"
fi

# =============================================================================
# SECTION 6: MetaEuk reference database — fungi RefSeq proteomes
# MetaEuk predicts eukaryotic gene models by matching contigs against a
# reference protein database. A comprehensive set of fungal RefSeq proteins
# covers ECM Basidiomycota (Amanita, Cortinarius, Suillus), Ascomycota
# saprotrophs (Trichoderma, Aspergillus), and Glomeromycota where available.
#
# Step 6a: Download all NCBI RefSeq fungi protein FASTAs via rsync
#          (fungi.1.protein.faa.gz, fungi.2.protein.faa.gz, ...)
# Step 6b: Concatenate and decompress into a single FASTA
# Step 6c: Build MMseqs2 database (MetaEuk uses MMseqs2 internally;
#          pre-building saves ~30 min per sample across 30-40 samples)
# =============================================================================

METAEUK_FASTA="${DB_DIR}/metaeuk/fungi_refseq_proteins.faa"
METAEUK_DB="${DB_DIR}/metaeuk/fungi_refseq_db"

if [[ ! -f "${METAEUK_FASTA}" ]]; then
    echo "--- [6a] Downloading NCBI RefSeq fungi proteins via rsync (~3-8 GB compressed) ---"
    # NCBI rsync host is rsync.ncbi.nlm.nih.gov (:: notation), NOT ftp.ncbi.nlm.nih.gov
    rsync -av --no-motd \
        --include="*.protein.faa.gz" \
        --exclude="*" \
        rsync.ncbi.nlm.nih.gov::refseq/release/fungi/ \
        "${SCRATCH_DIR}/fungi_refseq/"

    echo "--- [6b] Concatenating and decompressing fungi protein FASTAs ---"
    zcat "${SCRATCH_DIR}/fungi_refseq/"*.protein.faa.gz > "${METAEUK_FASTA}"
    echo "Fungi RefSeq FASTA done: $(date)"
else
    echo "[SKIP] Fungi RefSeq FASTA already exists"
fi

if [[ ! -f "${METAEUK_DB}" ]]; then
    echo "--- [6c] Building MMseqs2 database for MetaEuk ---"
    mmseqs createdb \
        "${METAEUK_FASTA}" \
        "${METAEUK_DB}" \
        --threads ${THREADS}
    echo "MetaEuk MMseqs2 database done: $(date)"
else
    echo "[SKIP] MetaEuk MMseqs2 database already exists"
fi

# =============================================================================
# Completion
# =============================================================================

echo ""
echo "============================================================"
echo "Database setup COMPLETE : $(date)"
echo ""
echo "Summary of database locations:"
echo "  MMseqs2 UniRef90 : ${DB_DIR}/mmseqs_taxonomy/uniref90"
echo "  NCBI taxonomy    : ${DB_DIR}/taxonomy/"
echo "  KOfam            : ${DB_DIR}/kofam/"
echo "  dbCAN3           : ${DB_DIR}/dbcan/"
echo "  PHI-base         : ${DB_DIR}/phibase/phi-base.dmnd"
echo "  MetaEuk target   : ${DB_DIR}/metaeuk/fungi_refseq_db"
echo "  MetaEuk FASTA    : ${DB_DIR}/metaeuk/fungi_refseq_proteins.faa"
echo "============================================================"

touch "${DB_DIR}/databases_complete.flag"
