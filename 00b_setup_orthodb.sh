#!/bin/bash
# =============================================================================
# 00b_setup_orthodb.sh
# Slurm job: build the OrthoDB Fungi near-single-copy ortholog DIAMOND
# database, used by 01b_orthodb_genecount.sh for genome-equivalent
# normalization (additive alternative to the Asparaginase/PF01112 method —
# see 08_integrate.R Section 6c and the README).
#
# Pulled out as a STANDALONE script (separate from 00_setup_databases.sh)
# specifically so this large, not-yet-fully-verified download/filter chain
# can be tested in isolation — a failure here (e.g. an unexpected FASTA
# header format) doesn't cost re-running the rest of database setup.
# Idempotent like every other script in this pipeline: re-running skips the
# build entirely if the final DIAMOND database already exists.
#
# Submit directly:
#   sbatch 00b_setup_orthodb.sh
# or run interactively on a compute node for faster iteration while
# debugging (the heavy steps are curl/zcat/awk/python, not GPU/MPI work).
#
# Uses the Fungi level (NCBI taxid 4751), NOT Dikarya (451864) — Dikarya is
# a valid NCBI taxon but OrthoDB does NOT compute orthologous groups at that
# specific level (confirmed directly: /search?level=451864 returns zero
# results even with no phyloprofile filter at all, while level=4751 returns
# real data). odb12v2_levels.tab lists exactly which levels OrthoDB has OGs
# for — check that file directly if retargeting this to a different clade.
# Fungi-level is arguably a better fit anyway: broader than Dikarya (also
# covers Glomeromycota/Chytridiomycota etc.), relevant since AM fungi may be
# present in non-ectomycorrhizal host samples (e.g. ash/elm) that Dikarya
# specifically would exclude.
#
# Uses OrthoDB v12v2's bulk tab-delimited data dump
# (https://data.orthodb.org/v12/download/odb_data_dump/), NOT the /search
# REST API — even with the taxid issue above resolved, the bulk dump is also
# what OrthoDB's own documentation recommends for this kind of large local
# processing, and avoids the API's 1 request/second rate limit entirely.
#
# The "near-single-copy, present in >90% of species" filter (matching the
# source paper's criterion and the web UI's Phyloprofile filter) is computed
# LOCALLY here from the dump's plain tab files, rather than relying on any
# server-side filter:
#   odb12v2_levels.tab.gz     (~16 KB)  -> total species count under Fungi
#   odb12v2_OGs.tab.gz        (~128 MB) -> OG id -> level (filter to Fungi)
#   odb12v2_OG2genes.tab.gz   (~4.5 GB) -> OG id -> gene id, filtered to
#                                          Fungi OGs only
#   odb12v2_genes.tab.gz      (~4.5 GB) -> gene id -> organism id, filtered
#                                          to genes from the step above
#   odb12v2_og_aa_fasta.gz    (~35 GB)  -> protein sequences, filtered down
#                                          to only the genes in kept OGs
#
# This is a large one-time download (~44 GB total) — comparable in scale to
# the UniRef90 build in 00_setup_databases.sh. Each large file is deleted
# immediately after being filtered down to free scratch space.
#
# NOTE — the FASTA header format in odb12v2_og_aa_fasta.gz was not confirmed
# ahead of implementation (OrthoDB's dump README only says headers contain
# "orthodb internal gene id as well as a public id"). Step [9f]'s matching
# logic is defensive (checks any whitespace/pipe-delimited token in the
# header) but should be spot-checked on first run — if N_SEQS comes back 0
# or suspiciously low, inspect a few header lines directly:
#   zcat <scratch>/orthodb_dump/og_aa_fasta.gz | head -50 | grep '^>'
# (note: og_aa_fasta.gz is deleted after this script runs successfully —
# re-download a small byte range with curl -r if you need to inspect it
# after a failure further down the pipeline.)
# =============================================================================

#SBATCH --job-name=metaG_setup_orthodb
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu
#SBATCH --output=logs/setup/setup_orthodb_%j.out
#SBATCH --error=logs/setup/setup_orthodb_%j.err

# NOTE: --output/--error above are relative to the directory you run `sbatch`
# from — that directory must already contain logs/setup/ before submitting:
#   mkdir -p logs/setup && sbatch 00b_setup_orthodb.sh
# (SLURM resolves --output/--error before the script body runs, so a mkdir
# inside this script would be too late to help the very first lines of output.)

set -euo pipefail

# ── Parameters (mirror 00_setup_databases.sh — keep these in sync) ───────────
THREADS=16
DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"
SCRATCH_DIR="/scratch.global/falb0011"

# ── Environment ───────────────────────────────────────────────────────────────
set +u
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_annotation
set -u

mkdir -p "${DB_DIR}/orthodb"

echo "============================================================"
echo "OrthoDB Fungi database setup started : $(date)"
echo "DB_DIR                               : ${DB_DIR}"
echo "SCRATCH_DIR                          : ${SCRATCH_DIR}"
echo "============================================================"

FUNGI_TAXID=4751
ORTHODB_DUMP_BASE="https://data.orthodb.org/v12/download/odb_data_dump"
ORTHODB_TMP="${SCRATCH_DIR}/orthodb_dump"
ORTHODB_FASTA="${DB_DIR}/orthodb/fungi_orthologs.faa"
ORTHODB_GENE2OG="${DB_DIR}/orthodb/gene2og.tsv"
ORTHODB_DB="${DB_DIR}/orthodb/fungi_orthologs"

if [[ -f "${ORTHODB_DB}.dmnd" ]]; then
    echo "[SKIP] OrthoDB Fungi database already exists: ${ORTHODB_DB}.dmnd"
    exit 0
fi

mkdir -p "${ORTHODB_TMP}"

echo "--- [1] Total Fungi species count (odb12v2_levels.tab) ---"
curl -s -o "${ORTHODB_TMP}/levels.tab.gz" "${ORTHODB_DUMP_BASE}/odb12v2_levels.tab.gz"
# levels.tab columns: level_taxid, name, gene_count, OG_count, species_count
TOTAL_FUNGI_SPECIES=$(zcat "${ORTHODB_TMP}/levels.tab.gz" \
    | awk -F'\t' -v t="${FUNGI_TAXID}" '$1==t {print $5}')
rm -f "${ORTHODB_TMP}/levels.tab.gz"
if [[ -z "${TOTAL_FUNGI_SPECIES}" ]]; then
    echo "ERROR: Fungi level (taxid ${FUNGI_TAXID}) not found in odb12v2_levels.tab" >&2
    exit 1
fi
echo "  Total Fungi species: ${TOTAL_FUNGI_SPECIES}"

echo "--- [2] Fungi-level OG ids (odb12v2_OGs.tab) ---"
curl -s -o "${ORTHODB_TMP}/OGs.tab.gz" "${ORTHODB_DUMP_BASE}/odb12v2_OGs.tab.gz"
# OGs.tab columns: OG_id, level_taxid, OG_name
zcat "${ORTHODB_TMP}/OGs.tab.gz" \
    | awk -F'\t' -v t="${FUNGI_TAXID}" '$2==t {print $1}' \
    > "${ORTHODB_TMP}/fungi_og_ids.txt"
rm -f "${ORTHODB_TMP}/OGs.tab.gz"
N_FUNGI_OGS=$(wc -l < "${ORTHODB_TMP}/fungi_og_ids.txt")
echo "  Fungi-level OGs found: ${N_FUNGI_OGS}"
if [[ "${N_FUNGI_OGS}" -lt 10 ]]; then
    echo "ERROR: suspiciously few Fungi OGs (${N_FUNGI_OGS}) — check" >&2
    echo "  FUNGI_TAXID and the odb12v2_OGs.tab column layout by hand." >&2
    exit 1
fi

echo "--- [3] OG -> gene membership for Fungi OGs (odb12v2_OG2genes.tab) ---"
curl -s -o "${ORTHODB_TMP}/OG2genes.tab.gz" "${ORTHODB_DUMP_BASE}/odb12v2_OG2genes.tab.gz"
# OG2genes.tab columns: OG_id, gene_id. Stream-filter to Fungi OGs only —
# avoids holding the full (multi-GB decompressed) file on disk.
zcat "${ORTHODB_TMP}/OG2genes.tab.gz" | awk -F'\t' -v idfile="${ORTHODB_TMP}/fungi_og_ids.txt" '
    BEGIN { while ((getline line < idfile) > 0) keep[line] = 1 }
    ($1 in keep) { print }
' > "${ORTHODB_TMP}/fungi_og2genes.tab"
rm -f "${ORTHODB_TMP}/OG2genes.tab.gz"
cut -f2 "${ORTHODB_TMP}/fungi_og2genes.tab" | sort -u > "${ORTHODB_TMP}/fungi_gene_ids.txt"
echo "  Fungi OG-gene pairs: $(wc -l < "${ORTHODB_TMP}/fungi_og2genes.tab")"
echo "  Distinct genes        : $(wc -l < "${ORTHODB_TMP}/fungi_gene_ids.txt")"

echo "--- [4] Gene -> organism mapping for those genes (odb12v2_genes.tab) ---"
curl -s -o "${ORTHODB_TMP}/genes.tab.gz" "${ORTHODB_DUMP_BASE}/odb12v2_genes.tab.gz"
# genes.tab columns: gene_id, organism_id, orig_seq_id, ... (10 cols total)
zcat "${ORTHODB_TMP}/genes.tab.gz" | awk -F'\t' -v idfile="${ORTHODB_TMP}/fungi_gene_ids.txt" '
    BEGIN { while ((getline line < idfile) > 0) keep[line] = 1 }
    ($1 in keep) { print $1 "\t" $2 }
' > "${ORTHODB_TMP}/fungi_gene2org.tab"
rm -f "${ORTHODB_TMP}/genes.tab.gz"

echo "--- [5] Computing near-single-copy filter (present >=90%, single-copy >=90%) ---"
# Replicates the source paper's Phyloprofile criterion locally: for each
# Fungi OG, count distinct organisms (presence) and how many of those
# organisms contribute exactly one gene (single-copy), rather than relying
# on OrthoDB's server-side filter (which wasn't usable — see header note).
python3 - "${ORTHODB_TMP}/fungi_og2genes.tab" "${ORTHODB_TMP}/fungi_gene2org.tab" \
         "${TOTAL_FUNGI_SPECIES}" "${ORTHODB_TMP}/kept_og2genes.tab" <<'PYEOF'
import sys
from collections import defaultdict

og2genes_file, gene2org_file, total_species_str, out_file = sys.argv[1:5]
total_species = float(total_species_str)

gene2org = {}
with open(gene2org_file) as f:
    for line in f:
        gene_id, org_id = line.rstrip("\n").split("\t")
        gene2org[gene_id] = org_id

og2genes = defaultdict(list)
with open(og2genes_file) as f:
    for line in f:
        og_id, gene_id = line.rstrip("\n").split("\t")
        og2genes[og_id].append(gene_id)

kept = 0
with open(out_file, "w") as out:
    for og_id, genes in og2genes.items():
        org_counts = defaultdict(int)
        for g in genes:
            org = gene2org.get(g)
            if org is not None:
                org_counts[org] += 1
        n_orgs = len(org_counts)
        if n_orgs == 0:
            continue
        n_single = sum(1 for c in org_counts.values() if c == 1)
        presence_frac = n_orgs / total_species
        singlecopy_frac = n_single / n_orgs
        if presence_frac >= 0.9 and singlecopy_frac >= 0.9:
            kept += 1
            for g in genes:
                out.write(og_id + "\t" + g + "\n")

sys.stderr.write("Kept %d near-single-copy Fungi OGs (out of %d)\n" % (kept, len(og2genes)))
PYEOF

N_KEPT_OGS=$(cut -f1 "${ORTHODB_TMP}/kept_og2genes.tab" | sort -u | wc -l)
echo "  Near-single-copy OGs kept: ${N_KEPT_OGS}"
if [[ "${N_KEPT_OGS}" -lt 50 ]]; then
    echo "WARNING: fewer than 50 OGs passed the near-single-copy filter." >&2
    echo "  Consider relaxing the 0.9/0.9 thresholds in this script's Python block if this seems too strict." >&2
fi
cut -f2 "${ORTHODB_TMP}/kept_og2genes.tab" | sort -u > "${ORTHODB_TMP}/kept_gene_ids.txt"

echo "--- [6] Extracting protein sequences for kept genes (odb12v2_og_aa_fasta, ~35 GB) ---"
# Idempotency: skip re-downloading if a previous (e.g. interrupted/failed) run
# already fetched this file — it's the single biggest, slowest download here.
if [[ -s "${ORTHODB_TMP}/og_aa_fasta.gz" ]]; then
    echo "  [SKIP] og_aa_fasta.gz already present in scratch, reusing it"
else
    curl -s -o "${ORTHODB_TMP}/og_aa_fasta.gz" "${ORTHODB_DUMP_BASE}/odb12v2_og_aa_fasta.gz"
fi

# Write the extraction script to a real file rather than a heredoc on the same
# line as the zcat pipe. `zcat f | python3 - <<'EOF' ... EOF` is broken: the
# heredoc and the pipe both target python3's stdin, and the heredoc wins (it's
# how `python3 -` gets its source code) — so the piped FASTA data never
# reaches the script's own `sys.stdin` at all, which is immediately at EOF by
# the time the for-loop runs. That produced "Wrote 0 matched sequences" with
# no error, despite the target IDs genuinely being present in the data
# (confirmed by direct grep during debugging). Writing to a separate .py file
# and invoking `python3 script.py` keeps stdin free for the actual pipe.
cat > "${ORTHODB_TMP}/extract_fasta.py" <<'PYEOF'
import sys, re

gene_ids_file, og2genes_file, fasta_out, gene2og_out = sys.argv[1:5]

target_genes = set()
with open(gene_ids_file) as f:
    for line in f:
        g = line.strip()
        if g:
            target_genes.add(g)

# OrthoDB gene id -> OG id (may be many-to-many if a gene is in multiple kept OGs)
gene_to_ogs = {}
with open(og2genes_file) as f:
    for line in f:
        og_id, gene_id = line.rstrip("\n").split("\t")
        gene_to_ogs.setdefault(gene_id, []).append(og_id)

n_written = 0
with open(fasta_out, "w") as fasta_f, open(gene2og_out, "w") as map_f:
    header = None
    seq_lines = []
    matched_id = None
    matched_ogs = None

    def flush():
        global n_written
        if header is not None and matched_id is not None:
            fasta_f.write(header + "\n")
            for sl in seq_lines:
                fasta_f.write(sl + "\n")
            for og in matched_ogs:
                map_f.write(matched_id + "\t" + og + "\n")
            n_written += 1

    for line in sys.stdin:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            header = line
            seq_lines = []
            tokens = re.split(r"[\s|]+", line[1:])
            matched_id = None
            matched_ogs = None
            for tok in tokens:
                if tok in target_genes:
                    matched_id = tokens[0]  # first token = what DIAMOND uses as sseqid
                    matched_ogs = gene_to_ogs.get(tok, [])
                    break
        else:
            if matched_id is not None:
                seq_lines.append(line)
    flush()

sys.stderr.write("Wrote %d matched sequences\n" % n_written)
PYEOF

zcat "${ORTHODB_TMP}/og_aa_fasta.gz" | python3 "${ORTHODB_TMP}/extract_fasta.py" \
    "${ORTHODB_TMP}/kept_gene_ids.txt" "${ORTHODB_TMP}/kept_og2genes.tab" \
    "${ORTHODB_FASTA}" "${ORTHODB_GENE2OG}"
rm -f "${ORTHODB_TMP}/og_aa_fasta.gz"

N_SEQS=$(grep -c '^>' "${ORTHODB_FASTA}" || echo 0)
echo "  Sequences extracted: ${N_SEQS}"
if [[ ! -s "${ORTHODB_FASTA}" ]]; then
    echo "ERROR: no sequences extracted — check FASTA header format assumptions" >&2
    echo "  in this script's Python block against a sample of the og_aa_fasta file." >&2
    exit 1
fi

echo "--- [7] Building DIAMOND database ---"
# DIAMOND not on PATH by default on this conda env's base shell on some
# nodes — load explicitly to be safe, matching 00_setup_databases.sh.
module load diamond/2.0.15-gcc-8.2.0-gkldzx7 2>/dev/null || true
diamond makedb --in "${ORTHODB_FASTA}" --db "${ORTHODB_DB}" --threads ${THREADS}

rm -rf "${ORTHODB_TMP}"

echo "============================================================"
echo "OrthoDB Fungi database COMPLETE : $(date)"
echo "  FASTA      : ${ORTHODB_FASTA}"
echo "  gene2og    : ${ORTHODB_GENE2OG}"
echo "  DIAMOND DB : ${ORTHODB_DB}.dmnd"
echo "============================================================"
