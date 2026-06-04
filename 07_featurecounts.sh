#!/bin/bash
# =============================================================================
# 07_featurecounts.sh
# Slurm array job: count reads per predicted gene using featureCounts.
#
# featureCounts (from the subread package) intersects each sample's sorted BAM
# (reads mapped to contigs in the assembly step) with MetaEuk's gene coordinate
# GFF to count how many reads overlap each predicted gene. This produces the
# raw count matrix used for differential abundance analysis in R.
#
# The BAM files from the assembly step mapped reads to all contigs. MetaEuk's
# GFF only has coordinates for eukaryotic genes. featureCounts counts reads at
# those gene positions and ignores reads on unannotated contigs — no re-mapping
# is needed.
#
# GFF feature type: MetaEuk v6 produces GFF3 with feature type "CDS" and
# gene model IDs in the "ID" attribute. featureCounts is configured to:
#   -t CDS    count reads overlapping CDS features
#   -g ID     group (summarize) counts by the ID attribute (one row per gene)
#
# NOTE: If featureCounts warns about missing or unmatched feature types, check
# the MetaEuk GFF with: grep -v '^#' ${SAMPLE}.gff | cut -f3 | sort -u
# and update -t below to match the feature type present in the file.
#
# Inputs:
#   ${ASSEMBLY_DIR}/coverage/bam/${SAMPLE}_sorted.bam   (from assembly step)
#   ${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.gff   (MetaEuk gene coords)
#
# Outputs:
#   ${ANNOTATION_DIR}/featurecounts/${SAMPLE}_counts.txt        raw count table
#   ${ANNOTATION_DIR}/featurecounts/${SAMPLE}_counts.txt.summary  QC summary
#
# Do not run directly. Submit via submit_featurecounts.sh.
# =============================================================================

#SBATCH --job-name=metaG_featurecounts
#SBATCH --partition=msismall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=falb0011@umn.edu

set -euo pipefail

# ── Conda (featureCounts / subread is in metaG_annotation environment) ────────
source /common/software/install/migrated/anaconda/python3-2020.07-mamba/etc/profile.d/conda.sh
conda activate metaG_annotation

# ── Parameters ────────────────────────────────────────────────────────────────
THREADS=8

# ── Resolve sample ────────────────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${ACTIVE_LIST:-${ANNOTATION_DIR}/sample_list.txt}")

BAM="${ASSEMBLY_DIR}/coverage/bam/${SAMPLE}_sorted.bam"
GFF="${ANNOTATION_DIR}/metaeuk/${SAMPLE}/${SAMPLE}.gff"
OUT_COUNTS="${ANNOTATION_DIR}/featurecounts/${SAMPLE}_counts.txt"

echo "============================================================"
echo "Sample     : ${SAMPLE}"
echo "Job ID     : ${SLURM_JOB_ID}  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start      : $(date)"
echo "============================================================"

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM file not found: ${BAM}" >&2
    echo "  Expected in the assembly output directory from step 4/5 of the assembly pipeline." >&2
    exit 1
fi
if [[ ! -f "${GFF}" ]]; then
    echo "ERROR: MetaEuk GFF not found: ${GFF}" >&2
    echo "  Run step 02 (submit_metaeuk.sh) first." >&2
    exit 1
fi

# If the GFF is empty (no eukaryotic genes predicted), write an empty count file
if [[ ! -s "${GFF}" ]] || ! grep -qv '^#' "${GFF}" 2>/dev/null; then
    echo "WARNING: GFF has no feature annotations for ${SAMPLE} — writing empty count file."
    mkdir -p "$(dirname "${OUT_COUNTS}")"
    touch "${OUT_COUNTS}"
    exit 0
fi

mkdir -p "$(dirname "${OUT_COUNTS}")"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: featureCounts
# -t CDS       : count reads overlapping CDS features (MetaEuk feature type)
# -g ID        : group counts by the gene model ID in the ID= attribute
# -O           : allow reads to be assigned to overlapping features (multiple
#                CDS features can be adjacent in MetaEuk output)
# --fraction   : if a read overlaps multiple features, distribute fractionally
#                (avoids double-counting where gene models overlap)
# -s 0         : unstranded counting — metagenome DNA reads have no strand info
# -p           : paired-end mode (our assembly used paired reads)
# --countReadPairs : count read pairs rather than individual reads
# -T           : number of threads
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "${OUT_COUNTS}" && -s "${OUT_COUNTS}" ]]; then
    echo "[SKIP] featureCounts output already exists"
else
    echo "--- Step 1: Detecting MetaEuk GFF structure ---"

    # MetaEuk v6 uses feature type "CDS" but older builds used "MetaEuk_CDS".
    # Auto-detect the most common non-comment feature type in the GFF.
    FEATURE_TYPE=$(grep -v '^#' "${GFF}" | cut -f3 | sort | uniq -c | sort -rn | awk 'NR==1{print $2}')
    echo "  GFF feature type detected  : ${FEATURE_TYPE}"

    # For multi-exon gene aggregation, prefer the Parent= attribute (links each
    # CDS/exon to its parent gene model). Fall back to ID= if Parent= is absent
    # (single-exon predictions where each CDS IS the gene model).
    if grep -v '^#' "${GFF}" | cut -f9 | grep -q 'Parent='; then
        GROUP_ATTR="Parent"
        echo "  Grouping attribute         : Parent (multi-exon hierarchy present)"
    else
        GROUP_ATTR="ID"
        echo "  Grouping attribute         : ID (no Parent attribute found)"
    fi

    echo "--- Step 1: featureCounts ---"
    featureCounts \
        -t "${FEATURE_TYPE}" \
        -g "${GROUP_ATTR}" \
        -O \
        --fraction \
        -s 0 \
        -p \
        --countReadPairs \
        -T ${THREADS} \
        -a "${GFF}" \
        -o "${OUT_COUNTS}" \
        "${BAM}"
    echo "Step 1 done: $(date)"

    # Sanity check: warn if almost no reads were assigned (likely wrong -t or -g)
    ASSIGNED=$(grep '^Assigned' "${OUT_COUNTS}.summary" | awk '{print $2}')
    TOTAL=$(awk 'NR>1 {sum+=$2} END {print sum+0}' "${OUT_COUNTS}.summary")
    if [[ -n "${ASSIGNED}" && -n "${TOTAL}" && "${TOTAL}" -gt 0 ]]; then
        PCT=$(awk "BEGIN{printf \"%.1f\", ${ASSIGNED}/${TOTAL}*100}")
        echo "  Read assignment rate: ${ASSIGNED}/${TOTAL} (${PCT}%)"
        if (( $(echo "${PCT} < 1.0" | awk '{print ($1 < 1.0)}') )); then
            echo "  WARNING: <1% reads assigned. Verify -t and -g match your GFF:"
            echo "    grep -v '^#' ${GFF} | cut -f3 | sort -u   # feature types"
            echo "    grep -v '^#' ${GFF} | cut -f9 | head -5   # attribute format"
        fi
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Summary
# The .summary file from featureCounts breaks down reads by assignment status:
# Assigned, Unassigned_NoFeatures, Unassigned_Ambiguity, etc.
# A low Assigned fraction is expected — most reads map to prokaryotic contigs
# or unannotated regions, not eukaryotic gene models.
# ─────────────────────────────────────────────────────────────────────────────

echo "--- Step 2: Assignment summary ---"

if [[ -f "${OUT_COUNTS}.summary" ]]; then
    echo "  featureCounts assignment summary for ${SAMPLE}:"
    cat "${OUT_COUNTS}.summary"
fi

N_GENES=$(grep -v '^#' "${OUT_COUNTS}" | grep -v '^Geneid' | wc -l || echo 0)
echo "  Gene models in count table : ${N_GENES}"

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "COMPLETE   : ${SAMPLE}"
echo "End        : $(date)"
echo "Key outputs:"
echo "  Count table : ${OUT_COUNTS}"
echo "  QC summary  : ${OUT_COUNTS}.summary"
echo "============================================================"
