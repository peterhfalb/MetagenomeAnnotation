#!/bin/bash
# =============================================================================
# run_annotation_pipeline.sh
# Master submission script: set your paths once, submit the full pipeline.
#
# Submits all steps in dependency order. Steps 3-7 run in parallel after
# MetaEuk (step 2) completes. Database setup (step 0) is skipped automatically
# if databases are already built.
#
# Completion detection: before each step is submitted, all samples are checked
# for existing output files. If every sample already has output, that step is
# skipped entirely (no SLURM job submitted). This makes re-running after a
# partial failure cheap — only the steps with missing output are resubmitted.
#
# Usage:
#   bash run_annotation_pipeline.sh              # normal run
#   bash run_annotation_pipeline.sh --test 2     # test on first 2 samples
#   bash run_annotation_pipeline.sh --force-db   # force re-run database setup
#
# If a step fails mid-run, see the RECOVERY section at the bottom of this file.
#
# After all jobs complete, run step 8 manually (Rscript 08_integrate.R).
# The exact command is printed at submission time.
# =============================================================================

set -euo pipefail

# =============================================================================
# !! EDIT THESE TWO PATHS BEFORE RUNNING !!
# =============================================================================

# Output directory from the MetagenomeAssembly pipeline
ASM_DIR="/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Assembled"

# Where annotation outputs will be written (created automatically)
ANN_DIR="/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation"

# =============================================================================
# Configuration — leave as-is unless you have a reason to change
# =============================================================================

DB_DIR="/projects/standard/kennedyp/shared/databases/metaG_annotation"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_N=""
FORCE_DB=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test)     TEST_N="${2:-2}"; shift 2 ;;
        --force-db) FORCE_DB=true;   shift   ;;
        *) echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

TEST_FLAG=""
[[ -n "${TEST_N}" ]] && TEST_FLAG="--test ${TEST_N}"

# ── Validate paths ────────────────────────────────────────────────────────────
if [[ "${ASM_DIR}" == "/path/to/MetagenomeAssembly/output" || \
      "${ANN_DIR}" == "/path/to/annotation/output" ]]; then
    echo "ERROR: Edit ASM_DIR and ANN_DIR at the top of this script before running." >&2
    exit 1
fi
if [[ ! -d "${ASM_DIR}" ]]; then
    echo "ERROR: Assembly directory not found: ${ASM_DIR}" >&2
    exit 1
fi

echo "============================================================"
echo "Annotation pipeline submission"
echo "Assembly dir    : ${ASM_DIR}"
echo "Annotation dir  : ${ANN_DIR}"
[[ -n "${TEST_N}" ]] && echo "TEST MODE       : first ${TEST_N} samples"
echo "============================================================"
echo ""

# =============================================================================
# STEP 0: Database setup — auto-skipped if already complete
# =============================================================================

DB_FLAG="${DB_DIR}/databases_complete.flag"
STEP0_ID=""

check_databases() {
    local all_present=true
    echo "  Database status:"

    # Parallel indexed arrays avoid the associative-array + set -u interaction
    # that causes "unbound variable" errors in bash 4.x.
    local -a labels=(
        "MMseqs2 UniRef90"
        "NCBI taxonomy"
        "KOfam profiles"
        "dbCAN3"
        "Pfam-A"
        "MetaEuk target DB"
    )
    local -a paths=(
        "${DB_DIR}/mmseqs_taxonomy/uniref90"
        "${DB_DIR}/taxonomy/nodes.dmp"
        "${DB_DIR}/kofam/profiles"
        "${DB_DIR}/dbcan/dbCAN.hmm.h3i"
        "${DB_DIR}/pfam/Pfam-A.hmm.h3i"
        "${DB_DIR}/metaeuk/fungi_refseq_db"
    )

    local i
    for i in "${!labels[@]}"; do
        if [[ -e "${paths[$i]}" ]]; then
            echo "    [OK]      ${labels[$i]}"
        else
            echo "    [MISSING] ${labels[$i]}"
            all_present=false
        fi
    done

    # PHI-base: FASTA requires manual download; DIAMOND DB is built by setup job.
    local phibase_fasta
    phibase_fasta=$(find "${DB_DIR}/phibase" -maxdepth 1 -name "*.fas" 2>/dev/null | head -1)
    if [[ -f "${DB_DIR}/phibase/phi-base.dmnd" ]]; then
        echo "    [OK]      PHI-base"
    elif [[ -n "${phibase_fasta}" ]]; then
        echo "    [OK]      PHI-base  (FASTA present — DIAMOND DB will be built by setup job)"
    else
        echo "    [MISSING] PHI-base  (FASTA not found — manual download required before setup job runs)"
        echo "              See Section 5 in 00_setup_databases.sh for instructions."
        all_present=false
    fi

    # Mycocosm/Phytozome: FASTA + taxid mapping require manual download; the
    # custom MMseqs2 taxonomy DB is built by the setup job. WARNING, not a
    # hard failure — this taxonomy layer is additive, the pipeline runs fine
    # without it (06_mmseqs_taxonomy.sh / UniRef90 still provides taxonomy).
    local myco_fasta myco_mapping
    myco_fasta=$(find "${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome" -maxdepth 1 -name "combined_proteins.faa" 2>/dev/null | head -1)
    myco_mapping="${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome/taxid_mapping.tsv"
    if [[ -f "${DB_DIR}/mmseqs_taxonomy/mycocosm_phytozome/db" ]]; then
        echo "    [OK]      Mycocosm/Phytozome taxonomy DB"
    elif [[ -n "${myco_fasta}" && -f "${myco_mapping}" ]]; then
        echo "    [OK]      Mycocosm/Phytozome taxonomy DB  (FASTA + mapping present — DB will be built by setup job)"
    else
        echo "    [MISSING] Mycocosm/Phytozome taxonomy DB  (manual download + taxid mapping required — optional, additive taxonomy layer)"
        echo "              See Section 8 in 00_setup_databases.sh for instructions."
        # Not added to all_present — this database is optional/additive, so its
        # absence should not block the rest of the pipeline from running.
    fi

    ${all_present}
}

if [[ "${FORCE_DB}" == true ]]; then
    echo "Step 0: --force-db set — submitting database setup regardless of current state."
    check_databases || true
    echo ""
    STEP0_ID=$(bash "${SCRIPT_DIR}/submit_setup_databases.sh" | grep -oP '(?<=job: )\d+')
    echo "  Step 0 job ID: ${STEP0_ID}"
    echo ""
elif [[ -f "${DB_FLAG}" ]]; then
    echo "Step 0: databases already complete (found ${DB_FLAG}) — skipping."
    echo ""
else
    echo "Step 0: databases_complete.flag not found — checking individual databases..."
    if check_databases; then
        echo ""
        echo "  All databases present but flag missing. Touching flag and skipping setup."
        touch "${DB_FLAG}"
    else
        echo ""
        echo "  One or more databases missing — submitting database setup job."
        echo "  NOTE: PHI-base requires a manual download first (see 00_setup_databases.sh Section 5)."
        echo "  If PHI-base is the only missing database, touch ${DB_FLAG} and re-run with --force-db=false."
        echo ""
        STEP0_ID=$(bash "${SCRIPT_DIR}/submit_setup_databases.sh" | grep -oP '(?<=job: )\d+')
        echo "  Step 0 job ID: ${STEP0_ID}"
    fi
    echo ""
fi

# =============================================================================
# Completion detection
# =============================================================================
# Before each step, check whether every sample already has its output file.
# If so, skip submission entirely (no SLURM job queued for that step).
# The check uses ${ASM_DIR}/sample_list.txt, which is always present after
# the assembly pipeline runs. For --test N, only the first N samples are checked.

CHECK_LIST=""
if [[ -f "${ASM_DIR}/sample_list.txt" ]]; then
    if [[ -n "${TEST_N}" ]]; then
        CHECK_LIST="$(mktemp)"
        head -n "${TEST_N}" "${ASM_DIR}/sample_list.txt" > "${CHECK_LIST}"
        trap 'rm -f "${CHECK_LIST}"' EXIT
    else
        CHECK_LIST="${ASM_DIR}/sample_list.txt"
    fi
fi

# step_complete STEP_NAME OUTPUT_PATTERN
# Returns 0 if all samples in CHECK_LIST have an existing output file.
# OUTPUT_PATTERN should contain the literal word SAMPLE (replaced per sample).
# Patterns with SAMPLE appearing more than once (e.g. metaeuk/SAMPLE/SAMPLE.fas)
# are handled correctly — all occurrences are substituted.
step_complete() {
    local step_name="$1"
    local pattern="$2"
    [[ -z "${CHECK_LIST}" ]] && return 1

    local n_done=0 n_missing=0
    while IFS= read -r s; do
        [[ -z "${s}" ]] && continue
        local out="${pattern//SAMPLE/${s}}"
        # compgen -G handles both exact paths and glob patterns (e.g. overview.*)
        if compgen -G "${out}" > /dev/null 2>&1; then
            n_done=$(( n_done + 1 ))
        else
            n_missing=$(( n_missing + 1 ))
        fi
    done < "${CHECK_LIST}"

    if [[ "${n_missing}" -eq 0 ]]; then
        echo "  [ALL DONE] ${step_name} — all ${n_done} sample(s) complete, skipping submission"
        return 0
    else
        echo "  [SUBMIT]   ${step_name} — ${n_missing}/$(( n_done + n_missing )) sample(s) need processing"
        return 1
    fi
}

# =============================================================================
# STEP 1: Tiara contig classification
# =============================================================================

echo "Step 1: Tiara contig classification"
STEP1_ID=""
if ! step_complete "Tiara" "${ANN_DIR}/tiara/SAMPLE_tiara.txt"; then
    STEP1_ARGS=(
        --assembly-dir   "${ASM_DIR}"
        --annotation-dir "${ANN_DIR}"
    )
    [[ -n "${TEST_FLAG}" ]] && STEP1_ARGS+=( ${TEST_FLAG} )
    [[ -n "${STEP0_ID}"  ]] && STEP1_ARGS+=( --after "${STEP0_ID}" )

    STEP1_ID=$(bash "${SCRIPT_DIR}/submit_classify_contigs.sh" "${STEP1_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP1_ID}"
fi
echo ""

# =============================================================================
# STEP 2: MetaEuk gene prediction
# =============================================================================

echo "Step 2: MetaEuk gene prediction"
STEP2_ID=""
if ! step_complete "MetaEuk" "${ANN_DIR}/metaeuk/SAMPLE/SAMPLE.fas"; then
    STEP2_ARGS=(
        --assembly-dir   "${ASM_DIR}"
        --annotation-dir "${ANN_DIR}"
    )
    [[ -n "${TEST_FLAG}" ]] && STEP2_ARGS+=( ${TEST_FLAG} )
    [[ -n "${STEP1_ID}"  ]] && STEP2_ARGS+=( --after "${STEP1_ID}" )

    STEP2_ID=$(bash "${SCRIPT_DIR}/submit_metaeuk.sh" "${STEP2_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP2_ID}"
fi
echo ""

# =============================================================================
# STEPS 3-7: Functional annotation + quantification (all parallel after step 2)
# =============================================================================

# Determine which upstream job the parallel steps should wait on.
# If step 2 was submitted, wait for it. If step 2 was skipped but step 1 was
# submitted, wait for step 1. If both were skipped, no dependency needed.
PARALLEL_DEP=""
if [[ -n "${STEP2_ID}" ]]; then
    PARALLEL_DEP="${STEP2_ID}"
elif [[ -n "${STEP1_ID}" ]]; then
    PARALLEL_DEP="${STEP1_ID}"
fi

echo "Steps 3-7: functional annotation + quantification (parallel)"
[[ -n "${PARALLEL_DEP}" ]] && echo "  Dependency: afterok:${PARALLEL_DEP}"
echo ""

PARALLEL_ARGS=( --annotation-dir "${ANN_DIR}" )
[[ -n "${TEST_FLAG}"    ]] && PARALLEL_ARGS+=( ${TEST_FLAG} )
[[ -n "${PARALLEL_DEP}" ]] && PARALLEL_ARGS+=( --after "${PARALLEL_DEP}" )

STEP3_ID=""
echo "Step 3: KOfamScan"
if ! step_complete "KOfamScan" "${ANN_DIR}/kofam/SAMPLE_kofam.tsv"; then
    STEP3_ID=$(bash "${SCRIPT_DIR}/submit_kofam.sh" "${PARALLEL_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP3_ID}"
fi
echo ""

STEP3B_ID=""
echo "Step 3b: Pfam domain annotation (HMMER)"
if ! step_complete "Pfam" "${ANN_DIR}/pfam/SAMPLE_pfam_mapper.tsv"; then
    STEP3B_ID=$(bash "${SCRIPT_DIR}/submit_pfam.sh" "${PARALLEL_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP3B_ID}"
fi
echo ""

STEP4_ID=""
echo "Step 4: dbCAN3"
if ! step_complete "dbCAN3" "${ANN_DIR}/dbcan/SAMPLE/overview.*"; then
    STEP4_ID=$(bash "${SCRIPT_DIR}/submit_dbcan.sh" "${PARALLEL_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP4_ID}"
fi
echo ""

STEP5_ID=""
echo "Step 5: PHI-base"
if ! step_complete "PHI-base" "${ANN_DIR}/phibase/SAMPLE_phibase.tsv"; then
    STEP5_ID=$(bash "${SCRIPT_DIR}/submit_phibase.sh" "${PARALLEL_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP5_ID}"
fi
echo ""

STEP6_ID=""
echo "Step 6: MMseqs2 taxonomy"
if ! step_complete "MMseqs2 taxonomy" "${ANN_DIR}/mmseqs_taxonomy/SAMPLE_lca.tsv"; then
    STEP6_ID=$(bash "${SCRIPT_DIR}/submit_mmseqs_taxonomy.sh" "${PARALLEL_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP6_ID}"
fi
echo ""

STEP6C_ID=""
echo "Step 6c: Mycocosm/Phytozome taxonomy (additive — no filtering applied)"
if ! step_complete "Mycocosm/Phytozome taxonomy" "${ANN_DIR}/mycocosm_taxonomy/SAMPLE_lca.tsv"; then
    STEP6C_ID=$(bash "${SCRIPT_DIR}/submit_mycocosm_taxonomy.sh" "${PARALLEL_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP6C_ID}"
fi
echo ""

STEP7_ID=""
echo "Step 7: featureCounts"
FC_ARGS=( --assembly-dir "${ASM_DIR}" --annotation-dir "${ANN_DIR}" )
[[ -n "${TEST_FLAG}"    ]] && FC_ARGS+=( ${TEST_FLAG} )
[[ -n "${PARALLEL_DEP}" ]] && FC_ARGS+=( --after "${PARALLEL_DEP}" )
if ! step_complete "featureCounts" "${ANN_DIR}/featurecounts/SAMPLE_counts.txt"; then
    STEP7_ID=$(bash "${SCRIPT_DIR}/submit_featurecounts.sh" "${FC_ARGS[@]}" \
        | grep -oP '(?<=job: )\d+')
    echo "  Job ID: ${STEP7_ID}"
fi
echo ""

# =============================================================================
# Summary
# =============================================================================

# Collect only the IDs of steps that were actually submitted.
SUBMITTED_IDS=()
[[ -n "${STEP0_ID}"  ]] && SUBMITTED_IDS+=("${STEP0_ID}")
[[ -n "${STEP1_ID}"  ]] && SUBMITTED_IDS+=("${STEP1_ID}")
[[ -n "${STEP2_ID}"  ]] && SUBMITTED_IDS+=("${STEP2_ID}")
[[ -n "${STEP3_ID}"  ]] && SUBMITTED_IDS+=("${STEP3_ID}")
[[ -n "${STEP3B_ID}" ]] && SUBMITTED_IDS+=("${STEP3B_ID}")
[[ -n "${STEP4_ID}"  ]] && SUBMITTED_IDS+=("${STEP4_ID}")
[[ -n "${STEP5_ID}"  ]] && SUBMITTED_IDS+=("${STEP5_ID}")
[[ -n "${STEP6_ID}"  ]] && SUBMITTED_IDS+=("${STEP6_ID}")
[[ -n "${STEP6C_ID}" ]] && SUBMITTED_IDS+=("${STEP6C_ID}")
[[ -n "${STEP7_ID}"  ]] && SUBMITTED_IDS+=("${STEP7_ID}")

if [[ "${#SUBMITTED_IDS[@]}" -eq 0 ]]; then
    echo "============================================================"
    echo "All steps already complete — no jobs submitted."
    echo ""
    echo "To force a step to re-run, delete its output files for the"
    echo "relevant sample(s) and rerun this script."
    echo ""
    echo "Step 8 — run now:"
    echo "  Rscript ${SCRIPT_DIR}/08_integrate.R \\"
    echo "    --annotation-dir ${ANN_DIR} \\"
    echo "    --assembly-dir   ${ASM_DIR} \\"
    echo "    --db-dir         ${DB_DIR}"
    echo "============================================================"
    exit 0
fi

ALL_IDS=$(IFS=','; echo "${SUBMITTED_IDS[*]}")

echo "============================================================"
echo "Jobs submitted:"
echo ""
echo "  Step 0 databases    : ${STEP0_ID:-skipped}"
echo "  Step 1 Tiara        : ${STEP1_ID:-skipped}"
echo "  Step 2 MetaEuk      : ${STEP2_ID:-skipped}"
echo "  Step 3 KOfamScan    : ${STEP3_ID:-skipped}"
echo "  Step 3b Pfam        : ${STEP3B_ID:-skipped}"
echo "  Step 4 dbCAN3       : ${STEP4_ID:-skipped}"
echo "  Step 5 PHI-base     : ${STEP5_ID:-skipped}"
echo "  Step 6 MMseqs2 tax  : ${STEP6_ID:-skipped}"
echo "  Step 6c Myco/Phyto  : ${STEP6C_ID:-skipped}"
echo "  Step 7 featureCounts: ${STEP7_ID:-skipped}"
echo ""
echo "Monitor:"
echo "  squeue -j ${ALL_IDS}"
echo "  tail -f ${ANN_DIR}/logs/..."
echo ""
echo "Logs: ${ANN_DIR}/logs/"
echo ""
echo "Step 8 — run after all jobs complete:"
echo "  Rscript ${SCRIPT_DIR}/08_integrate.R \\"
echo "    --annotation-dir ${ANN_DIR} \\"
echo "    --assembly-dir   ${ASM_DIR} \\"
echo "    --db-dir         ${DB_DIR}"
echo "============================================================"
echo ""
# =============================================================================
# RECOVERY — what to do if a step fails
# =============================================================================
# SLURM dependency type is 'afterok', meaning: if any task in an array job
# exits non-zero, all downstream dependent jobs are cancelled with status
# 'DependencyNeverSatisfied'. The completion checks in this script mean that
# re-running after a partial failure is safe — completed samples are detected
# and skipped, so only the failed samples are resubmitted.
#
# To recover after a failure:
#
# 1. Find failed samples:
#      grep -l "FAILED\|Error\|error" ${ANN_DIR}/logs/<step>/*.err
#      or check which output files are missing:
#      for s in $(cat ${ANN_DIR}/sample_list.txt); do
#          [[ ! -f ${ANN_DIR}/metaeuk/${s}/${s}.fas ]] && echo "MISSING: ${s}"
#      done
#
# 2. Simply rerun this script — it will detect completed samples and only
#    submit jobs for the steps/samples that are still missing output.
#
# 3. Alternatively, rerun only specific samples for a specific step:
#      bash submit_metaeuk.sh \
#        --assembly-dir ${ASM_DIR} --annotation-dir ${ANN_DIR} \
#        --samples "FAILED_SAMPLE_01,FAILED_SAMPLE_02"
#
# Note: steps 3-7 are independent of each other. If only step 6 (MMseqs2 taxonomy)
# fails, steps 3-5 and 7 are unaffected and you only need to rerun step 6.
#
# Steps 3b (Pfam) and 6c (Mycocosm/Phytozome taxonomy) are likewise independent
# of the rest of the parallel block and of each other — both depend only on
# step 2 (MetaEuk). Step 6c is additive/optional: if its database isn't built
# (see Section 8 in 00_setup_databases.sh), check_databases() warns but does
# not block the rest of the pipeline from running.
