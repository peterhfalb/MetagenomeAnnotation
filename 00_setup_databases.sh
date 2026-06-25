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
#   pfam/         Pfam-A HMM profiles (domain annotation via hmmscan)
#   mmseqs_taxonomy/mycocosm_phytozome/  Custom MMseqs2 taxonomy DB built from
#                 user-supplied Mycocosm (fungal) + Phytozome (plant) proteins
#
# NOTE — PHI-base requires a manual download step before this job runs.
#        See SECTION 5 below for instructions.
#
# NOTE — Mycocosm/Phytozome requires a manual download + taxid mapping step
#        before this job runs. See SECTION 8 below for instructions.
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
set +u
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_annotation
set -u

# DIAMOND is still needed for the PHI-base and fungi DIAMOND database builds.
# BLAST+ module no longer required (replaced by mmseqs databases for taxonomy).
module load diamond/2.0.15-gcc-8.2.0-gkldzx7

# ── Create directories ────────────────────────────────────────────────────────
mkdir -p "${DB_DIR}"/{taxonomy,mmseqs_taxonomy,diamond,kofam,dbcan,phibase,metaeuk,pfam}
mkdir -p "${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome"
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
# SECTION 4: dbCAN — CAZyme databases
# Installed tool is the standalone `dbcan` package v5.2.9, a ground-up
# rewrite of the older v3/v4 run_dbcan — different CLI, different DB layout.
# run_dbcan operates in three modes (04_dbcan.sh --methods hmm,diamond,dbCANsub):
#   hmm       — dbCAN HMM profiles (most sensitive, primary result)
#   diamond   — searched against characterized CAZy proteins (used as evidence)
#   dbCANsub  — substrate-annotated CAZy subfamily HMMs (glycan substrate call)
# Genes are called high-confidence CAZymes when hmm and diamond agree — this
# version doesn't compute that agreement itself, so 04_dbcan.sh/08_integrate.R
# do it by comparing gene IDs across the two result files.
# `run_dbcan database` downloads everything needed, including the CGC/PUL
# bacterial gene-cluster files (TF.hmm, STP.hmm, dbCAN-PUL.xlsx, etc.) that
# this pipeline doesn't use — confirmed to download successfully on this
# version/site without the 403 issues seen on older dbCAN releases, so they
# are no longer stubbed; the stub fallback below is left in only as a
# defensive no-op in case that changes.
# =============================================================================

# NOTE: the substrate HMM file is named "dbCAN-sub.hmm" (hyphen), not
# "dbCAN_sub.hmm" (underscore) — easy to mis-type; double-check both sides of
# any future edits here against the actual filename in ${DB_DIR}/dbcan/.
DBCAN_SUB_HMM="${DB_DIR}/dbcan/dbCAN-sub.hmm"

if [[ ! -f "${DB_DIR}/dbcan/dbCAN.hmm.h3i" ]]; then
    echo "--- [4] Downloading dbCAN databases ---"
    run_dbcan database --db_dir "${DB_DIR}/dbcan" || true

    # Verify the files we actually need downloaded successfully
    if [[ ! -f "${DB_DIR}/dbcan/CAZy.dmnd" ]]; then
        echo "ERROR: CAZy.dmnd not found after dbCAN download — check network access." >&2
        exit 1
    fi
    if [[ ! -f "${DB_DIR}/dbcan/dbCAN.hmm" ]]; then
        echo "ERROR: dbCAN.hmm not found after dbCAN download — check network access." >&2
        exit 1
    fi
    if [[ ! -f "${DBCAN_SUB_HMM}" ]]; then
        echo "ERROR: ${DBCAN_SUB_HMM} not found after dbCAN download — check network access." >&2
        exit 1
    fi

    # Press both HMM databases (may not have run if download errored after
    # file transfer, or — for dbCAN-sub specifically — because run_dbcan
    # database doesn't press it itself, only downloads the raw file).
    hmmpress "${DB_DIR}/dbcan/dbCAN.hmm"
    hmmpress "${DBCAN_SUB_HMM}"

    # Defensive fallback only — confirmed unnecessary on this site/version
    # (see comment above), kept in case a future environment 403s again.
    for opt_file in TF.hmm STP.hmm dbCAN-PUL.xlsx; do
        [[ -f "${DB_DIR}/dbcan/${opt_file}" ]] || touch "${DB_DIR}/dbcan/${opt_file}"
    done

    if [[ ! -s "${DB_DIR}/dbcan/fam-substrate-mapping.tsv" ]]; then
        echo "WARNING: fam-substrate-mapping.tsv missing or empty after run_dbcan"
        echo "  database — dbCAN-sub substrate names may come back blank or as raw"
        echo "  subfamily IDs instead of substrate names. Download it manually and"
        echo "  place it at ${DB_DIR}/dbcan/fam-substrate-mapping.tsv"
    fi

    echo "dbCAN done: $(date)"
else
    echo "[SKIP] dbCAN databases already exist"
fi

# dbCAN-sub.hmm may not have been pressed in a previous run of this script
# (added after the initial dbCAN setup, and run_dbcan database does not press
# it itself) — check independently of the main dbCAN.hmm.h3i sentinel above
# so re-running this script picks it up.
if [[ -f "${DBCAN_SUB_HMM}" && ! -f "${DBCAN_SUB_HMM}.h3i" ]]; then
    echo "--- [4b] Pressing dbCAN-sub.hmm (substrate prediction) ---"
    hmmpress "${DBCAN_SUB_HMM}"
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

PHIBASE_DMND="${DB_DIR}/phibase/phi-base.dmnd"

# Accept any .fas file in the phibase directory (filename varies by release)
PHIBASE_FASTA=$(find "${DB_DIR}/phibase" -maxdepth 1 -name "*.fas" | head -1)

if [[ -z "${PHIBASE_FASTA}" ]]; then
    echo ""
    echo "WARNING: PHI-base FASTA not found — skipping PHI-base DIAMOND build."
    echo "  Download the current FASTA from https://www.phi-base.org/ and place in:"
    echo "  ${DB_DIR}/phibase/"
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
    echo "--- [6a] Downloading NCBI RefSeq fungi proteins via HTTPS (~3-8 GB compressed) ---"
    # rsync port 873 is blocked on MSI compute nodes; wget -r is also unreliable because
    # set -e kills the script if wget returns non-zero on any link error. Instead:
    # step 1 — fetch directory listing; step 2 — download each file individually.
    mkdir -p "${SCRATCH_DIR}/fungi_refseq"

    echo "  Fetching NCBI fungi RefSeq directory listing..."
    wget -O "${SCRATCH_DIR}/fungi_index.html" \
        "https://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/"

    grep -oP 'fungi\.\d+\.protein\.faa\.gz' "${SCRATCH_DIR}/fungi_index.html" \
        | sort -u > "${SCRATCH_DIR}/fungi_files.txt"

    N_EXPECTED=$(wc -l < "${SCRATCH_DIR}/fungi_files.txt")
    echo "  Found ${N_EXPECTED} protein FASTA files to download"
    if [[ "${N_EXPECTED}" -eq 0 ]]; then
        echo "ERROR: Could not parse fungi file list from NCBI — check network access." >&2
        exit 1
    fi

    while IFS= read -r fname; do
        dest="${SCRATCH_DIR}/fungi_refseq/${fname}"
        if [[ -f "${dest}" ]]; then
            echo "  [skip] ${fname}"
        else
            echo "  Downloading ${fname}..."
            wget -c -q "https://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/${fname}" \
                -O "${dest}"
        fi
    done < "${SCRATCH_DIR}/fungi_files.txt"

    N_FILES=$(ls "${SCRATCH_DIR}/fungi_refseq/"*.protein.faa.gz 2>/dev/null | wc -l)
    echo "  ${N_FILES} / ${N_EXPECTED} files present after download"
    if [[ "${N_FILES}" -eq 0 ]]; then
        echo "ERROR: No fungi protein files downloaded." >&2
        exit 1
    fi

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
# SECTION 7: Pfam-A — domain annotation HMM profiles
# Used by hmmscan (step 03b) to assign Pfam domains to predicted proteins via
# the family-specific "gathering threshold" (--cut_ga), Pfam's own curated
# cutoff rather than a single fixed E-value across all families. Also fetches
# Pfam-A.hmm.dat (family accession → human-readable name/description) so
# downstream tables aren't limited to bare PF##### accessions.
# =============================================================================

PFAM_HMM="${DB_DIR}/pfam/Pfam-A.hmm"

if [[ ! -f "${PFAM_HMM}.h3i" ]]; then
    echo "--- [7] Downloading and indexing Pfam-A (~100 MB compressed) ---"
    # If this 403s from a compute node (EBI FTP can be flaky from MSI compute
    # nodes — dbCAN's Section 4 hits the same class of issue), re-run just this
    # wget on the login node, then resubmit the job; hmmpress will pick up from
    # the already-downloaded file.
    wget -q -O "${PFAM_HMM}.gz" \
        https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip "${PFAM_HMM}.gz"

    wget -q -O "${DB_DIR}/pfam/Pfam-A.hmm.dat.gz" \
        https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
    gunzip "${DB_DIR}/pfam/Pfam-A.hmm.dat.gz"

    hmmpress "${PFAM_HMM}"
    echo "Pfam-A done: $(date)"
else
    echo "[SKIP] Pfam-A database already exists"
fi

# =============================================================================
# SECTION 8: Mycocosm + Phytozome — custom MMseqs2 taxonomy database
# Provides finer-resolution fungal taxonomy (down to subphylum, e.g.
# Agaricomycotina / Pezizomycotina) and flags plant-derived sequences, by
# searching against JGI's curated fungal (Mycocosm) and plant (Phytozome)
# proteomes. Runs ADDITIONALLY alongside the UniRef90 MMseqs2 taxonomy step
# (Section 2 above) — this does not replace it, and no filtering is applied
# anywhere in the pipeline based on this taxonomy; it only adds resolution.
#
# *** MANUAL STEPS REQUIRED BEFORE THIS JOB RUNS ***
# JGI genomes require a signed Data Use Agreement and cannot be auto-downloaded:
#   1. Register at https://genome.jgi.doe.gov/portal/ and accept the Mycocosm
#      and Phytozome data use policies.
#   2. Download protein FASTA files for your genome panel — recommended: fungal
#      genomes spanning Agaricomycotina and Pezizomycotina (the two EM-forming
#      subphyla), plus a panel of plant genomes for host/contamination flagging.
#   3. Concatenate all downloaded FASTAs into one file and place it at:
#        ${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome/combined_proteins.faa
#   4. Build a two-column TSV (no header) mapping each FASTA header ID (or a
#      unique prefix shared by all headers from one genome) to its NCBI taxid,
#      and place it at:
#        ${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome/taxid_mapping.tsv
#      This step is required because JGI headers do not carry NCBI taxonomy
#      IDs the way UniProt/UniRef headers do (unlike Section 2's UniRef90
#      build, which has NCBI taxonomy mapping built in automatically). Use
#      NCBI-recognized taxids only — taxonomizr rank expansion in step 08
#      depends on it.
# This script will then build the database automatically.
# =============================================================================

MYCO_DIR="${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome"
MYCO_FASTA=$(find "${MYCO_DIR}" -maxdepth 1 -name "combined_proteins.faa" | head -1)
MYCO_MAPPING="${MYCO_DIR}/taxid_mapping.tsv"
MYCO_DB="${MYCO_DIR}/db"

if [[ -z "${MYCO_FASTA}" || ! -f "${MYCO_MAPPING}" ]]; then
    echo ""
    echo "WARNING: Mycocosm/Phytozome FASTA or taxid mapping not found — skipping build."
    echo "  Required files:"
    echo "    ${MYCO_DIR}/combined_proteins.faa"
    echo "    ${MYCO_MAPPING}"
    echo "  See SECTION 8 comments above for download/mapping instructions."
    echo "  Then re-run this script; the database build will run automatically."
    echo ""
elif [[ ! -f "${MYCO_DB}" ]]; then
    echo "--- [8a] Validating taxid mapping coverage ---"
    # Catch unmapped headers now rather than letting createtaxdb silently
    # assign taxid 0 to anything missing from the mapping file.
    UNMAPPED=$(grep '^>' "${MYCO_FASTA}" | sed 's/^>//' | awk '{print $1}' \
        | sort -u > "${SCRATCH_DIR}/myco_headers.txt"; \
        cut -f1 "${MYCO_MAPPING}" | sort -u > "${SCRATCH_DIR}/myco_mapped.txt"; \
        comm -23 "${SCRATCH_DIR}/myco_headers.txt" "${SCRATCH_DIR}/myco_mapped.txt")

    if [[ -n "${UNMAPPED}" ]]; then
        echo "ERROR: Some FASTA header IDs have no entry in taxid_mapping.tsv:" >&2
        echo "${UNMAPPED}" | head -20 >&2
        echo "  (showing up to 20; see ${SCRATCH_DIR}/myco_headers.txt vs ${SCRATCH_DIR}/myco_mapped.txt)" >&2
        echo "  Fix taxid_mapping.tsv to cover all headers, then re-run." >&2
        exit 1
    fi
    echo "  All FASTA headers have a taxid mapping."

    echo "--- [8b] Building Mycocosm/Phytozome MMseqs2 taxonomy database ---"
    mmseqs createdb \
        "${MYCO_FASTA}" \
        "${MYCO_DB}" \
        --threads ${THREADS}

    # NOTE: verify exact flag names against the installed MMseqs2 version
    # (mmseqs createtaxdb -h) before relying on this in production — flag
    # naming for custom tax-mapping-file-driven taxonomy DBs has changed
    # across MMseqs2 releases.
    mmseqs createtaxdb \
        "${MYCO_DB}" \
        "${SCRATCH_DIR}/myco_taxdb_tmp" \
        --tax-mapping-file "${MYCO_MAPPING}" \
        --ncbi-tax-dump "${DB_DIR}/taxonomy"

    echo "Mycocosm/Phytozome taxonomy DB done: $(date)"
else
    echo "[SKIP] Mycocosm/Phytozome taxonomy database already exists"
fi

# =============================================================================
# Completion
# =============================================================================

echo ""
echo "============================================================"
echo "Database setup COMPLETE : $(date)"
echo ""
echo "Summary of database locations:"
echo "  MMseqs2 UniRef90       : ${DB_DIR}/mmseqs_taxonomy/uniref90"
echo "  Mycocosm/Phytozome     : ${MYCO_DB} $([[ -f "${MYCO_DB}" ]] && echo '[OK]' || echo '[MISSING - see SECTION 8]')"
echo "  NCBI taxonomy          : ${DB_DIR}/taxonomy/"
echo "  KOfam                  : ${DB_DIR}/kofam/"
echo "  dbCAN3                 : ${DB_DIR}/dbcan/"
echo "  PHI-base               : ${DB_DIR}/phibase/phi-base.dmnd"
echo "  Pfam-A                 : ${PFAM_HMM}"
echo "  MetaEuk target         : ${DB_DIR}/metaeuk/fungi_refseq_db"
echo "  MetaEuk FASTA          : ${DB_DIR}/metaeuk/fungi_refseq_proteins.faa"
echo "============================================================"

touch "${DB_DIR}/databases_complete.flag"
