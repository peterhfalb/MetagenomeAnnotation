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
SCRATCH_DIR="/users/9/falb0011/FAB2/Summer2022/HMSC/bacteria"

# BLAST+ module name for blastdbcmd (used to extract NR sequences from MSI copy)
# Run: module avail blast    — find the correct module name, then update below
BLAST_MODULE="blast-plus/2.13.0-gcc-8.2.0-vo4mr4d"

# ── Environment ───────────────────────────────────────────────────────────────
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_annotation

module load diamond/2.0.15-gcc-8.2.0-gkldzx7
module load "${BLAST_MODULE}"

# ── Create directories ────────────────────────────────────────────────────────
mkdir -p "${DB_DIR}"/{taxonomy,diamond,kofam,dbcan,phibase,metaeuk}
mkdir -p "${SCRATCH_DIR}"/{fungi_refseq}

echo "============================================================"
echo "Database setup started : $(date)"
echo "DB_DIR                 : ${DB_DIR}"
echo "SCRATCH_DIR            : ${SCRATCH_DIR}"
echo "============================================================"

# =============================================================================
# SECTION 1: NCBI Taxonomy
# Provides the species/lineage lookup tables used by DIAMOND (--taxonnodes,
# --taxonnames) and downstream R analysis for interpreting taxon IDs.
# prot.accession2taxid maps every NR protein accession to an NCBI taxon ID —
# required at DIAMOND database build time so per-protein taxonomy is embedded.
# =============================================================================

if [[ ! -f "${DB_DIR}/taxonomy/nodes.dmp" ]]; then
    echo "--- [1a] Downloading NCBI taxdump ---"
    wget -q -O "${SCRATCH_DIR}/taxdump.tar.gz" \
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf "${SCRATCH_DIR}/taxdump.tar.gz" -C "${DB_DIR}/taxonomy/"
    rm "${SCRATCH_DIR}/taxdump.tar.gz"
    echo "taxdump done: $(date)"
else
    echo "[SKIP] taxonomy/nodes.dmp already exists"
fi

if [[ ! -f "${DB_DIR}/taxonomy/prot.accession2taxid.gz" ]]; then
    echo "--- [1b] Downloading prot.accession2taxid (~5 GB compressed) ---"
    wget -q -O "${DB_DIR}/taxonomy/prot.accession2taxid.gz" \
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
    echo "prot.accession2taxid done: $(date)"
else
    echo "[SKIP] prot.accession2taxid already exists"
fi

# =============================================================================
# SECTION 2: DIAMOND NR database
# Built from MSI's existing quarterly-updated BLAST NR at
# /common/bioref/blast/latest/nr — avoids downloading ~400 GB of FASTA.
# blastdbcmd streams all sequences to stdout; diamond makedb reads from stdin.
# Taxonomy is embedded at build time using files from Section 1.
# Final .dmnd is ~80 GB and stored permanently.
# Runtime: ~10-14 hours on 32 threads.
# =============================================================================

NR_DMND="${DB_DIR}/diamond/nr.dmnd"

if [[ ! -f "${NR_DMND}" ]]; then
    echo "--- [2] Building DIAMOND NR database from MSI BLAST NR (~10-14 hr) ---"

    blastdbcmd \
        -db /common/bioref/blast/latest/nr \
        -entry all \
        -out /dev/stdout \
    | diamond makedb \
        --in /dev/stdin \
        --db "${NR_DMND}" \
        --taxonmap "${DB_DIR}/taxonomy/prot.accession2taxid.gz" \
        --taxonnodes "${DB_DIR}/taxonomy/nodes.dmp" \
        --taxonnames "${DB_DIR}/taxonomy/names.dmp" \
        --threads ${THREADS}

    echo "DIAMOND NR done: $(date)"
else
    echo "[SKIP] DIAMOND NR already exists (${NR_DMND})"
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
    cd "${DB_DIR}/dbcan"
    run_dbcan download_db "${DB_DIR}/dbcan"
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
echo "  DIAMOND NR       : ${DB_DIR}/diamond/nr.dmnd"
echo "  NCBI taxonomy    : ${DB_DIR}/taxonomy/"
echo "  KOfam            : ${DB_DIR}/kofam/"
echo "  dbCAN3           : ${DB_DIR}/dbcan/"
echo "  PHI-base         : ${DB_DIR}/phibase/phi-base.dmnd"
echo "  MetaEuk target   : ${DB_DIR}/metaeuk/fungi_refseq_db"
echo "  MetaEuk FASTA    : ${DB_DIR}/metaeuk/fungi_refseq_proteins.faa"
echo "============================================================"

touch "${DB_DIR}/databases_complete.flag"
