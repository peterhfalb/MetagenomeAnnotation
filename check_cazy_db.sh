#!/bin/bash
# =============================================================================
# check_cazy_db.sh
# Read-only summary of how sequences in the CAZy readmap database are labeled,
# to check what taxonomic information each header format carries. Uses the
# source FASTA if still on scratch, otherwise reads headers from the .dmnd.
#
# Outputs (in the output directory, default ~/cazy_db_check):
#   headers.txt        all sequence headers
#   summary.txt        header formats, examples, and taxonmap coverage
#   jgi_genomes.tsv    one row per JGI source file: code, source file,
#                      n_sequences, n_distinct_numbers. The number in JGI headers
#                      varies per sequence, so it can't identify an organism.
#   taxid_lineages.tsv one row per taxonmap taxid: taxid, group (from its full
#                      NCBI lineage), phylum, n_accessions
#
# Fungal status is taken from the lineage: a taxid is fungal only if its parent
# chain reaches Fungi (NCBI taxid 4751).
#
# Usage:
#   module load diamond/2.0.15-gcc-8.2.0-gkldzx7
#   bash check_cazy_db.sh [output_dir]
# =============================================================================

set -euo pipefail

DB_DIR="${DB_DIR:-/projects/standard/kennedyp/shared/databases/metaG_annotation/cazy_readmap}"
SOURCE_FASTA="${SOURCE_FASTA:-/scratch.global/falb0011/CAZyDB.fa}"
OUT_DIR="${1:-${HOME}/cazy_db_check}"
TAXONMAP="${DB_DIR}/cazy_taxonmap.tsv"
NODES="${NODES:-${DB_DIR}/nodes_patched.dmp}"
NAMES="${NAMES:-${DB_DIR}/../taxonomy/names.dmp}"
PARSER_PHYLA="Ascomycota,Basidiomycota,Chytridiomycota,Mucoromycota,Glomeromycota,Blastocladiomycota,Neocallimastigomycota,Zoopagomycota,Mortierellomycota,Microsporidia"
LABEL_RE='^(GH|GT|PL|CE|AA|CBM)[0-9]+(_[0-9]+)?$'

mkdir -p "${OUT_DIR}"
HEADERS="${OUT_DIR}/headers.txt"

# ── Collect headers ───────────────────────────────────────────────────────────
if [[ -s "${SOURCE_FASTA}" ]]; then
    echo "Reading headers from ${SOURCE_FASTA}"
    LC_ALL=C grep -a '^>' "${SOURCE_FASTA}" | sed 's/^>//' > "${HEADERS}"
else
    echo "Source FASTA not found; reading headers from the .dmnd"
    if ! command -v diamond > /dev/null; then
        echo "ERROR: diamond not on PATH. Load it first, then re-run:" >&2
        echo "  module load diamond/2.0.15-gcc-8.2.0-gkldzx7" >&2
        exit 1
    fi
    diamond getseq --db "${DB_DIR}/CAZy_taxonomy" --quiet \
        | LC_ALL=C grep -a '^>' | sed 's/^>//' > "${HEADERS}"
fi
if [[ ! -s "${HEADERS}" ]]; then
    echo "ERROR: no headers collected in ${HEADERS}" >&2
    exit 1
fi

# ── Summarize ─────────────────────────────────────────────────────────────────
{
    echo "=== Total sequences ==="
    wc -l < "${HEADERS}"

    echo; echo "=== Sequences by header format ==="
    awk -F'|' -v re="${LABEL_RE}" '
        { f = ($1 ~ re) ? "JGI (family first)" : (($1 ~ /\./) ? "accession, versioned" : "accession, unversioned")
          n[f]++ }
        END { for (k in n) print n[k] "\t" k }' "${HEADERS}" | sort -nr

    echo; echo "=== Examples of each format ==="
    awk -F'|' -v re="${LABEL_RE}" '
        { f = ($1 ~ re) ? "JGI" : (($1 ~ /\./) ? "versioned" : "unversioned")
          if (seen[f]++ < 6) print f "\t" $0 }' "${HEADERS}"

    echo; echo "=== Non-family-first headers that are not accession|labels ==="
    awk -F'|' -v re="${LABEL_RE}" '
        $1 !~ re { for (i = 2; i <= NF; i++) if ($i !~ re && $i !~ /^[0-9]+\.[0-9-]+/) { print; break } }' \
        "${HEADERS}" | head -20

    echo; echo "=== JGI headers: last field shape (digits replaced, top 15) ==="
    awk -F'|' -v re="${LABEL_RE}" '
        $1 ~ re { x = $NF; sub(/^[^_]*_/, "", x); gsub(/[0-9]+/, "N", x); n[x]++ }
        END { for (k in n) print n[k] "\t" k }' "${HEADERS}" | sort -nr | head -15

    if [[ -s "${TAXONMAP}" ]]; then
        echo; echo "=== Taxonmap coverage of accession headers ==="
        awk -F'|' -v re="${LABEL_RE}" '
            NR == FNR { split($0, t, "\t"); if (FNR > 1) map[t[1]] = 1; next }
            $1 !~ re { total++; if ($1 in map) hit++ }
            END { printf "%d of %d accession headers have a taxid (%.2f%%)\n", hit, total, 100 * hit / total }' \
            "${TAXONMAP}" "${HEADERS}"
    fi
} > "${OUT_DIR}/summary.txt"

# ── JGI genomes ───────────────────────────────────────────────────────────────
{
    printf 'code\tsource_file\tn_sequences\tn_distinct_numbers\n'
    awk -F'|' -v re="${LABEL_RE}" '
        $1 ~ re {
            num = "NA"
            for (i = 2; i < NF; i++) if ($i ~ /^[0-9]+$/) num = $i
            n[$NF]++
            if (!(($NF, num) in seen)) { seen[$NF, num] = 1; d[$NF]++ }
        }
        END {
            for (f in n) { code = f; sub(/_.*$/, "", code); print code "\t" f "\t" n[f] "\t" d[f] }
        }' "${HEADERS}" | sort -t$'\t' -k1,1
} > "${OUT_DIR}/jgi_genomes.tsv"

# ── NCBI lineage of every taxonmap taxid ──────────────────────────────────────
if [[ -s "${TAXONMAP}" && -s "${NODES}" && -s "${NAMES}" ]]; then
    {
        printf 'taxid\tgroup\tphylum\tn_accessions\n'
        awk -F'\t' '
            FILENAME == ARGV[1] { parent[$1] = $3; rank[$1] = $5; next }
            FILENAME == ARGV[2] { if ($7 == "scientific name") name[$1] = $3; next }
            FNR > 1 { n[$2]++ }
            END {
                for (t in n) {
                    g = "other_or_unresolved"; phy = "NA"; x = t; steps = 0
                    if (!(t in parent)) g = "taxid_not_in_nodes"
                    while ((x in parent) && steps++ < 100) {
                        if (rank[x] == "phylum") phy = name[x]
                        if      (x == 4751)  { g = "Fungi" }
                        else if (x == 2)     { g = "Bacteria" }
                        else if (x == 2157)  { g = "Archaea" }
                        else if (x == 10239) { g = "Viruses" }
                        else if (x == 2759 && g != "Fungi") { g = "Eukaryota_nonfungal" }
                        if (x == parent[x]) break
                        x = parent[x]
                    }
                    print t "\t" g "\t" phy "\t" n[t]
                }
            }' "${NODES}" "${NAMES}" "${TAXONMAP}" | sort -t$'\t' -k2,2 -k3,3
    } > "${OUT_DIR}/taxid_lineages.tsv"

    {
        echo; echo "=== Taxonmap accessions by lineage group ==="
        awk -F'\t' 'NR > 1 { n[$2] += $4 } END { for (k in n) print n[k] "\t" k }' \
            "${OUT_DIR}/taxid_lineages.tsv" | sort -nr

        echo; echo "=== Fungal accessions by phylum (* = missing from FUNGAL_PHYLA in 08_integrate.R) ==="
        awk -F'\t' -v known="${PARSER_PHYLA}" '
            BEGIN { split(known, a, ","); for (i in a) k[a[i]] = 1 }
            NR > 1 && $2 == "Fungi" { n[$3] += $4 }
            END { for (p in n) print n[p] "\t" p ((p in k) ? "" : " *") }' \
            "${OUT_DIR}/taxid_lineages.tsv" | sort -nr
    } >> "${OUT_DIR}/summary.txt"
else
    echo "[SKIP] taxonmap, nodes, or names not found; lineage check not run." >> "${OUT_DIR}/summary.txt"
fi

cat "${OUT_DIR}/summary.txt"
echo
echo "JGI genome rows : $(($(wc -l < "${OUT_DIR}/jgi_genomes.tsv") - 1))"
echo "Distinct codes  : $(tail -n +2 "${OUT_DIR}/jgi_genomes.tsv" | cut -f1 | sort -u | wc -l)"
echo "Written to      : ${OUT_DIR}"
