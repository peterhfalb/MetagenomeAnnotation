#!/usr/bin/env Rscript
# =============================================================================
# 08_integrate.R
# Merge all per-sample annotation outputs into analysis-ready tables.
#
# Usage (command line):
#   Rscript 08_integrate.R \
#     --annotation-dir /path/to/annotation \
#     --assembly-dir   /path/to/assembly \
#     --db-dir         /projects/standard/kennedyp/shared/databases/metaG_annotation \
#     [--out-dir       /path/to/output]   defaults to annotation-dir/integrated
#     [--min-tools     2]                 dbCAN: minimum tools for high-confidence call
#     [--top-hits      10]                DIAMOND NR: top N hits per gene for LCA
#     [--bitscore-frac 0.90]              LCA: fraction of top bitscore to retain hits
#     [--pseudocount   1]                 ALR normalization: pseudocount added to
#                                          numerator/denominator (see Section 6b)
#
# Interactive use (in an R session):
#   Set variables in the CONFIG block below, then source() this file.
#
# Outputs (all in out-dir/):
#   gene_annotations.tsv               one row per gene — contig coords + all annotation
#                                       layers (KO, CAZy, PHI, Pfam, UniRef90 + Mycocosm/
#                                       Phytozome taxonomy)
#   gene_counts_raw.tsv                raw count matrix: genes × samples (sparse — one
#                                       assembly per sample means each gene appears in
#                                       one sample only)
#   ko_count_matrix.tsv                KO × sample count matrix       ← DESeq2 / vegan
#   cazy_count_matrix.tsv              CAZy family × sample matrix    ← DESeq2 / vegan
#   phi_count_matrix.tsv               PHI-base phenotype × sample    ← pathogenicity
#   pfam_count_matrix.tsv              Pfam family × sample count matrix
#   *_count_matrix_normalized.tsv      ALR-normalized versions of all four matrices
#                                       above (Asparaginase/PF01112 pseudo-genome
#                                       normalization — see Section 6b and README)
#   asparaginase_pseudogenome_counts.tsv  raw PF01112 read count per sample (the ALR
#                                       denominator) — audit/QC file
#   summary_stats.tsv                  per-sample QC: gene counts, annotation rates
#
# NOTE on count normalization:
#   The raw *_count_matrix.tsv outputs are RAW read counts. Normalize in R
#   before differential analysis:
#   - DESeq2: pass raw counts directly (DESeq2 handles its own normalization)
#   - vegan:  normalize with decostand() (e.g., "total", "log", "hellinger")
#   - Do NOT pre-normalize before DESeq2.
#   The *_normalized.tsv outputs are a SEPARATE, already-log-transformed
#   pseudo-genome (ALR) normalization — do not feed these into DESeq2/vegan,
#   which expect raw counts; they are for direct compositional comparison.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(optparse)
  library(taxonomizr)
})

# =============================================================================
# ARGUMENT PARSING
# (When running interactively, set these variables and skip the parse_args block)
# =============================================================================

option_list <- list(
  make_option("--annotation-dir", type = "character", default = NULL,
              help = "Annotation pipeline output directory (required)"),
  make_option("--assembly-dir",   type = "character", default = NULL,
              help = "Assembly pipeline output directory (required)"),
  make_option("--db-dir",         type = "character",
              default = "/projects/standard/kennedyp/shared/databases/metaG_annotation",
              help = "Path to databases dir [default: %default]"),
  make_option("--out-dir",        type = "character", default = NULL,
              help = "Output directory [default: annotation-dir/integrated]"),
  make_option("--min-tools",      type = "integer",   default = 2L,
              help = "dbCAN min tools for high-confidence CAZyme [default: %default]"),
  make_option("--top-hits",       type = "integer",   default = 10L,
              help = "DIAMOND NR: top N hits per gene for LCA [default: %default]"),
  make_option("--bitscore-frac",  type = "double",    default = 0.90,
              help = "LCA: retain hits within this fraction of top bitscore [default: %default]"),
  make_option("--pseudocount",    type = "double",    default = 1.0,
              help = "ALR normalization: pseudocount for numerator/denominator [default: %default]"),
  make_option("--metadata",       type = "character", default = NULL,
              help = "Path to stakes.csv with stake_id and stand_type columns (optional; merges stand_type into summary_stats.tsv)"),
  make_option("--cazy-pident",   type = "double",    default = 50.0,
              help = "CAZy readmap: minimum % identity for a hit to be counted [default: %default] (Bahram 2018)"),
  make_option("--cazy-evalue",   type = "double",    default = 1e-9,
              help = "CAZy readmap: maximum e-value for a hit to be counted [default: %default] (Bahram 2018)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt[["annotation-dir"]]) || is.null(opt[["assembly-dir"]])) {
  stop("--annotation-dir and --assembly-dir are required.", call. = FALSE)
}

ann_dir    <- opt[["annotation-dir"]]
asm_dir    <- opt[["assembly-dir"]]
db_dir     <- opt[["db-dir"]]
out_dir    <- if (is.null(opt[["out-dir"]])) file.path(ann_dir, "integrated") else opt[["out-dir"]]
min_tools  <- opt[["min-tools"]]
top_hits   <- opt[["top-hits"]]
bs_frac    <- opt[["bitscore-frac"]]
pseudocount   <- opt[["pseudocount"]]
metadata_file  <- opt[["metadata"]]
cazy_pident    <- opt[["cazy-pident"]]
cazy_evalue    <- opt[["cazy-evalue"]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Integration started :", format(Sys.time()), "\n")
cat("Annotation dir      :", ann_dir, "\n")
cat("Output dir          :", out_dir, "\n")
cat("============================================================\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Read featureCounts output for one sample.
# Returns data.table with columns: gene_id, contig, start, end, strand, length, count
read_featurecounts <- function(file, sample_name) {
  if (!file.exists(file) || file.info(file)$size == 0) return(NULL)
  # featureCounts prepends several comment lines starting with '#' before the
  # column header. data.table's skip = "#" would START reading from the first '#'
  # line, not skip past them. Use grep -v to strip comments before fread sees them.
  dt <- fread(cmd = paste0("grep -v '^#' ", shQuote(file)),
              header = TRUE, sep = "\t", quote = "")
  if (nrow(dt) == 0) return(NULL)
  # Column names: Geneid Chr Start End Strand Length <bam_path>
  count_col <- ncol(dt)
  result <- dt[, .(
    gene_id = normalize_gene_id(sub("_mRNA$", "", Geneid)),
    contig  = Chr,
    start   = Start,
    end     = End,
    strand  = Strand,
    length  = Length,
    count   = .SD[[count_col]],
    sample  = sample_name
  )]
  result
}

# Read KOfamScan detail-tsv output.
# Returns data.table: gene_id, KO, KO_definition (significant hits only, * in col 1)
read_kofam <- function(file) {
  empty <- data.table(gene_id = character(), KO = character(), KO_definition = character())
  if (!file.exists(file) || file.info(file)$size == 0) return(empty)
  dt <- fread(file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
  if (nrow(dt) == 0) return(empty)
  if (any(dt$V1 == "*")) {
    # exec_annotation detail-tsv: threshold_flag | gene_id | KO | threshold | score | evalue | definition
    sig <- dt[V1 == "*"]
    sig[, .(gene_id = normalize_gene_id(V2), KO = V3,
            KO_definition = if (ncol(dt) >= 7L) V7 else NA_character_)]
  } else {
    # mapper format: gene_id | KO  (no threshold flag column)
    dt[, .(gene_id = normalize_gene_id(V1), KO = V2, KO_definition = NA_character_)]
  }
}

# Read dbCAN output for one sample.
# Installed tool is the standalone `dbcan` package v5.2.9 — a ground-up
# rewrite of the older v3/v4 run_dbcan that this function originally
# targeted. It does NOT write a unified overview file; each evidence stream
# (--methods hmm,diamond,dbCANsub in 04_dbcan.sh) writes its own file, and
# there is no built-in "#ofTools" agreement count. This function reads all
# three files directly and computes the HMMER∩DIAMOND agreement itself.
#
# Returns data.table: gene_id, CAZy_HMMER, CAZy_DIAMOND, CAZy_n_tools,
#                      CAZy_substrate
# CAZy_n_tools counts only HMMER/DIAMOND agreement (the two family-calling
# tools) — dbCAN_sub's substrate call is reported, not voted on.
read_dbcan <- function(dbcan_dir) {
  empty <- data.table(gene_id = character(), CAZy_HMMER = character(),
                      CAZy_DIAMOND = character(), CAZy_n_tools = integer(),
                      CAZy_substrate = character())
  if (is.null(dbcan_dir) || !dir.exists(dbcan_dir)) return(empty)

  read_one <- function(f) {
    if (!file.exists(f) || file.info(f)$size == 0) return(NULL)
    dt <- fread(f, header = TRUE, sep = "\t", quote = "", fill = TRUE)
    if (nrow(dt) == 0) return(NULL)
    dt
  }

  # HMMER family calls: "HMM Name" (e.g. "AA6.hmm") per "Target Name" (gene).
  # A gene can have multiple domain hits; collapse to one row per gene,
  # joining distinct families with ";" (same multi-hit convention used for
  # KO/Pfam elsewhere in this script).
  hmm_dt <- read_one(file.path(dbcan_dir, "dbCAN_hmm_results.tsv"))
  hmm_collapsed <- if (!is.null(hmm_dt) &&
                        all(c("Target Name", "HMM Name") %in% names(hmm_dt))) {
    setnames(hmm_dt, c("Target Name", "HMM Name"), c("gene_id", "family"))
    hmm_dt[, .(gene_id = normalize_gene_id(gene_id),
               family  = sub("\\.hmm$", "", family))][
      , .(CAZy_HMMER = paste(unique(family), collapse = ";")), by = gene_id]
  } else data.table(gene_id = character(), CAZy_HMMER = character())

  # DIAMOND family calls: "CAZy ID" embeds family as "ACCESSION|FAMILY".
  diamond_dt <- read_one(file.path(dbcan_dir, "diamond.out"))
  diamond_collapsed <- if (!is.null(diamond_dt) &&
                            all(c("Gene ID", "CAZy ID") %in% names(diamond_dt))) {
    setnames(diamond_dt, c("Gene ID", "CAZy ID"), c("gene_id", "family"))
    diamond_dt[, .(gene_id = normalize_gene_id(gene_id),
                   family  = sub("^[^|]*\\|", "", family))][
      , .(CAZy_DIAMOND = paste(unique(family), collapse = ";")), by = gene_id]
  } else data.table(gene_id = character(), CAZy_DIAMOND = character())

  # dbCAN-sub substrate calls: "Substrate" column, already translated from
  # raw subfamily hits via fam-substrate-mapping.tsv internally by run_dbcan.
  sub_dt <- read_one(file.path(dbcan_dir, "dbCANsub_hmm_results.tsv"))
  sub_collapsed <- if (!is.null(sub_dt) &&
                        all(c("Target Name", "Substrate") %in% names(sub_dt))) {
    setnames(sub_dt, "Target Name", "gene_id")
    sub_dt[, .(gene_id = normalize_gene_id(gene_id), substrate = Substrate)][
      !is.na(substrate) & !substrate %in% c("", "-", "N/A"),
      .(CAZy_substrate = paste(unique(substrate), collapse = ";")), by = gene_id]
  } else data.table(gene_id = character(), CAZy_substrate = character())

  all_ids <- unique(c(hmm_collapsed$gene_id, diamond_collapsed$gene_id, sub_collapsed$gene_id))
  if (length(all_ids) == 0) return(empty)

  result <- data.table(gene_id = all_ids)
  result <- hmm_collapsed[result, on = "gene_id"]
  result <- diamond_collapsed[result, on = "gene_id"]
  result <- sub_collapsed[result, on = "gene_id"]

  result[, CAZy_n_tools := as.integer(!is.na(CAZy_HMMER)) + as.integer(!is.na(CAZy_DIAMOND))]
  result
}

# Parse PHI-base stitle into components.
# PHI-base headers follow: PHI:ACCESSION|gene|pathogen|host|phenotype|...
# or older format: gene_name pathogen_species | phenotype
# Returns a list with: phi_accession, phi_gene, phi_pathogen, phi_phenotype
parse_phi_stitle <- function(stitle) {
  # PHI-base stitle format (# delimited):
  #   UniProtID # PHI:XXXX # gene_name # taxid # pathogen # phenotype
  # Multiple PHI IDs / phenotypes are joined with __ within a field.
  # Take the first phenotype when multiple are present.
  if (grepl("#", stitle, fixed = TRUE)) {
    parts <- str_split(stitle, "#")[[1]]
    list(
      phi_accession = str_extract(parts[2], "PHI:[0-9]+"),
      phi_gene      = if (length(parts) >= 3) trimws(parts[3]) else NA_character_,
      phi_pathogen  = if (length(parts) >= 5) trimws(parts[5]) else NA_character_,
      phi_phenotype = if (length(parts) >= 6) trimws(strsplit(parts[6], "__", fixed = TRUE)[[1]][1]) else NA_character_
    )
  } else {
    list(phi_accession = str_extract(stitle, "PHI:[0-9]+"),
         phi_gene = NA_character_, phi_pathogen = NA_character_, phi_phenotype = NA_character_)
  }
}

# Read PHI-base DIAMOND hits.
# Returns data.table: gene_id, phi_accession, phi_gene, phi_pathogen, phi_phenotype,
#                     phi_pident, phi_evalue, phi_bitscore
read_phibase <- function(file) {
  empty <- data.table(gene_id = character(), phi_accession = character(),
                      phi_gene = character(), phi_pathogen = character(),
                      phi_phenotype = character(), phi_pident = numeric(),
                      phi_evalue = numeric(), phi_bitscore = numeric())
  if (!file.exists(file) || file.info(file)$size == 0) return(empty)
  dt <- fread(file, header = FALSE, sep = "\t", quote = "",
    col.names = c("gene_id","sseqid","pident","length","mismatch","gapopen",
                  "qstart","qend","sstart","send","evalue","bitscore","stitle"))
  if (nrow(dt) == 0) return(empty)
  # Parse stitle into structured columns
  parsed <- rbindlist(lapply(dt$stitle, parse_phi_stitle))
  cbind(dt[, .(gene_id = normalize_gene_id(gene_id), phi_pident = pident, phi_evalue = evalue, phi_bitscore = bitscore)],
        parsed)
}

# Read MMseqs2 easy-taxonomy _lca.tsv output.
# Returns data.table: gene_id, lca_taxid, lca_rank, lca_name, lineage
# MMseqs2 computes LCA internally — one row per gene, already resolved.
# lca_taxid 0 = unclassified (no hits above threshold), stored as NA.
read_mmseqs_taxonomy <- function(file) {
  empty <- data.table(gene_id = character(), lca_taxid = integer(),
                      lca_rank = character(), lca_name = character(),
                      lineage = character())
  if (!file.exists(file) || file.info(file)$size == 0) return(empty)
  dt <- fread(file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
  if (nrow(dt) == 0) return(empty)
  result <- data.table(
    gene_id   = normalize_gene_id(as.character(dt[[1]])),
    lca_taxid = as.integer(dt[[2]]),
    lca_rank  = if (ncol(dt) >= 3) as.character(dt[[3]]) else NA_character_,
    lca_name  = if (ncol(dt) >= 4) as.character(dt[[4]]) else NA_character_,
    lineage   = if (ncol(dt) >= 5) as.character(dt[[5]]) else NA_character_
  )
  result[lca_taxid == 0L, lca_taxid := NA_integer_]
  result
}

# Read MMseqs2 easy-taxonomy _lca.tsv output from the Mycocosm/Phytozome
# custom database (06c_mycocosm_taxonomy.sh). Identical schema to
# read_mmseqs_taxonomy() (UniRef90) — same function body, renamed so the
# distinction at the call site is unambiguous. Column names get a myco_
# prefix downstream (Section 5) to avoid colliding with the UniRef90 columns
# when both taxonomy layers are joined into the same master table.
read_mycocosm_taxonomy <- function(file) {
  empty <- data.table(gene_id = character(), lca_taxid = integer(),
                      lca_rank = character(), lca_name = character(),
                      lineage = character())
  if (!file.exists(file) || file.info(file)$size == 0) return(empty)
  dt <- fread(file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
  if (nrow(dt) == 0) return(empty)
  result <- data.table(
    gene_id   = normalize_gene_id(as.character(dt[[1]])),
    lca_taxid = as.integer(dt[[2]]),
    lca_rank  = if (ncol(dt) >= 3) as.character(dt[[3]]) else NA_character_,
    lca_name  = if (ncol(dt) >= 4) as.character(dt[[4]]) else NA_character_,
    lineage   = if (ncol(dt) >= 5) as.character(dt[[5]]) else NA_character_
  )
  result[lca_taxid == 0L, lca_taxid := NA_integer_]
  result
}

# Read Pfam hmmscan mapper output (03b_pfam.sh).
# Returns data.table: gene_id, Pfam_accession, Pfam_name (one row per domain
# hit — a gene with multiple distinct Pfam domains has multiple rows). Every
# row already passed that family's gathering threshold (--cut_ga), so unlike
# read_kofam() there is no significance flag to filter on.
read_pfam <- function(file) {
  empty <- data.table(gene_id = character(), Pfam_accession = character(),
                      Pfam_name = character())
  if (!file.exists(file) || file.info(file)$size == 0) return(empty)
  dt <- fread(file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
  if (nrow(dt) == 0) return(empty)
  dt[, .(gene_id = normalize_gene_id(V1), Pfam_name = V2, Pfam_accession = V3)]
}

# Normalize MetaEuk gene IDs to a consistent 4-field join key:
#   targetID | contigID | strand | lowerBound
#
# GFF3 IDs (featureCounts, 4 fields):  targetID|contig|strand|lowerBound
# FASTA IDs (annotation tools, 9 fields): targetID|contig|strand|score|evalue|numExons|lowerBound|...
#
# The two sources differ in which field holds lowerBound: field 4 in GFF3 vs
# field 7 in FASTA. Detect format by pipe count and extract accordingly.
normalize_gene_id <- function(ids) {
  n_pipes <- nchar(ids) - nchar(gsub("|", "", ids, fixed = TRUE))
  result  <- ids

  # FASTA format (6+ pipes = 7+ fields): lowerBound is at field 7
  fas <- n_pipes >= 5L
  if (any(fas)) {
    result[fas] <- sub(
      "^([^|]+)\\|([^|]+)\\|([^|]+)\\|[^|]+\\|[^|]+\\|[^|]+\\|([^|]+)(\\|.*)?$",
      "\\1|\\2|\\3|\\4",
      ids[fas]
    )
  }

  # GFF format (≤4 pipes = ≤5 fields): lowerBound is already at field 4
  gff <- !fas
  if (any(gff)) {
    result[gff] <- sub("^([^|]+\\|[^|]+\\|[^|]+\\|[^|]+)(\\|.*)?$", "\\1", ids[gff])
  }

  result
}

# =============================================================================
# SECTION 1: Load sample list
# =============================================================================

sample_list_file <- file.path(ann_dir, "sample_list.txt")
if (!file.exists(sample_list_file)) stop("Sample list not found: ", sample_list_file)
samples <- readLines(sample_list_file)
samples <- samples[nzchar(samples)]
cat("Samples to integrate:", length(samples), "\n\n")

# =============================================================================
# SECTION 2: Read all per-sample data
# =============================================================================

cat("--- Section 2: Reading per-sample annotation files ---\n")

counts_list  <- list()  # featureCounts per sample
kofam_list   <- list()  # KOfamScan per sample
dbcan_list   <- list()  # dbCAN3 per sample
phi_list     <- list()  # PHI-base per sample
mmseqs_list  <- list()  # MMseqs2 taxonomy (UniRef90) per sample
pfam_list    <- list()  # Pfam (hmmscan) per sample
myco_list    <- list()  # MMseqs2 taxonomy (Mycocosm/Phytozome) per sample

for (s in samples) {
  cat("  Reading:", s, "\n")

  # featureCounts
  counts_list[[s]] <- read_featurecounts(
    file.path(ann_dir, "featurecounts", paste0(s, "_counts.txt")), s)

  # KOfamScan (reads the full detail-tsv; mapper is derived from it in read_kofam)
  kofam_list[[s]] <- read_kofam(file.path(ann_dir, "kofam", paste0(s, "_kofam.tsv")))

  # dbCAN3 — find actual overview file (name varies across v4 sub-versions)
  dbcan_list[[s]] <- read_dbcan(file.path(ann_dir, "dbcan", s))

  # PHI-base
  phi_list[[s]] <- read_phibase(file.path(ann_dir, "phibase", paste0(s, "_phibase.tsv")))

  # MMseqs2 taxonomy vs UniRef90 (_lca.tsv — LCA already computed by MMseqs2)
  mmseqs_list[[s]] <- read_mmseqs_taxonomy(
    file.path(ann_dir, "mmseqs_taxonomy", paste0(s, "_lca.tsv")))

  # Pfam domain annotation (gene_id -> Pfam_accession mapper, --cut_ga hits)
  pfam_list[[s]] <- read_pfam(file.path(ann_dir, "pfam", paste0(s, "_pfam_mapper.tsv")))

  # MMseqs2 taxonomy vs Mycocosm/Phytozome (additive — finer fungal subphylum
  # resolution + plant-sequence flagging; absent if the optional DB wasn't built)
  myco_list[[s]] <- read_mycocosm_taxonomy(
    file.path(ann_dir, "mycocosm_taxonomy", paste0(s, "_lca.tsv")))
}

# Combine across samples
counts_all <- rbindlist(counts_list,  use.names = TRUE, fill = TRUE)
kofam_all  <- rbindlist(kofam_list,   use.names = TRUE, fill = TRUE)
dbcan_all  <- rbindlist(dbcan_list,   use.names = TRUE, fill = TRUE)
pfam_all   <- rbindlist(pfam_list,    use.names = TRUE, fill = TRUE)
myco_all   <- rbindlist(myco_list,    use.names = TRUE, fill = TRUE)
phi_all    <- rbindlist(phi_list,     use.names = TRUE, fill = TRUE)
mmseqs_all <- rbindlist(mmseqs_list,  use.names = TRUE, fill = TRUE)

# Rename myco_all's columns immediately so they never collide with
# mmseqs_all's identically-named columns (lca_taxid, lca_rank, lca_name,
# lineage) once both taxonomy layers are joined into base_dt (Section 5).
if (nrow(myco_all) > 0) {
  setnames(myco_all,
           c("lca_taxid", "lca_rank", "lca_name", "lineage"),
           c("myco_lca_taxid", "myco_lca_rank", "myco_lca_name", "myco_lineage"))
} else {
  myco_all <- data.table(gene_id = character(), myco_lca_taxid = integer(),
                          myco_lca_rank = character(), myco_lca_name = character(),
                          myco_lineage = character())
}

cat("\n  Genes with featureCounts data :", nrow(counts_all), "\n")
cat("  KOfam significant hits        :", nrow(kofam_all), "\n")
cat("  dbCAN gene calls              :", nrow(dbcan_all), "\n")
cat("  PHI-base hits                 :", nrow(phi_all), "\n")
cat("  MMseqs2 taxonomy assignments  :", sum(!is.na(mmseqs_all$lca_taxid)), "\n")
cat("  Pfam domain hits               :", nrow(pfam_all), "\n")
cat("  Mycocosm/Phytozome assignments :", sum(!is.na(myco_all$myco_lca_taxid)), "\n\n")

# =============================================================================
# SECTION 3: Build raw count matrix (gene × sample)
# =============================================================================

cat("--- Section 3: Building raw count matrix ---\n")

# NOTE: Because each gene exists in exactly one sample's assembly, this matrix
# is essentially sparse. The analysis-useful matrices are the family-level
# count matrices built in Section 6.

if (nrow(counts_all) > 0) {
  count_matrix <- dcast(
    counts_all[, .(gene_id, sample, count)],
    gene_id ~ sample,
    value.var = "count",
    fill = 0L
  )
  cat("  Gene count matrix dimensions:", nrow(count_matrix), "genes ×",
      ncol(count_matrix) - 1, "samples\n\n")
} else {
  warning("No featureCounts data found. Count matrix will be empty.")
  count_matrix <- data.table(gene_id = character())
}

# =============================================================================
# SECTION 4: Expand MMseqs2 LCA taxon IDs to standard ranks
# =============================================================================
# MMseqs2 already computed the LCA internally — one taxid per gene, resolved.
# This section only needs to expand each unique LCA taxid to standard ranks
# (superkingdom, phylum, class, order, family, genus, species) using taxonomizr.
# This is far simpler than the previous approach: one lookup per unique taxid
# rather than processing top-N hits per gene and computing LCA in R.
#
# Two taxonomy sources are expanded here, sharing one getTaxonomy() lookup to
# avoid redundant queries: the existing UniRef90 LCA (mmseqs_all) and the
# additive Mycocosm/Phytozome LCA (myco_all). The latter only resolves
# correctly if taxid_mapping.tsv (00_setup_databases.sh Section 8) used valid
# NCBI taxids — see README for this dependency.

cat("--- Section 4: Expanding LCA taxon IDs to standard ranks ---\n")

tax_sql   <- file.path(db_dir, "taxonomy", "taxonomy.sql")
names_dmp <- file.path(db_dir, "taxonomy", "names.dmp")
nodes_dmp <- file.path(db_dir, "taxonomy", "nodes.dmp")
tax_ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

lca_table      <- NULL  # UniRef90, columns prefixed tax_*
myco_lca_table <- NULL  # Mycocosm/Phytozome, columns prefixed myco_tax_*

all_taxids_present <- (nrow(mmseqs_all) > 0 && !all(is.na(mmseqs_all$lca_taxid))) ||
                       (nrow(myco_all)   > 0 && !all(is.na(myco_all$myco_lca_taxid)))

if (!file.exists(names_dmp) || !file.exists(nodes_dmp)) {
  warning("NCBI taxonomy files not found — skipping rank expansion.\n",
          "  Expected: ", names_dmp)
} else if (!all_taxids_present) {
  warning("No classified MMseqs2 hits (UniRef90 or Mycocosm/Phytozome) — skipping rank expansion.")
} else {

  # Build SQLite taxonomy DB from local names.dmp/nodes.dmp (once, ~5-10 min)
  if (!file.exists(tax_sql)) {
    cat("  Building taxonomizr SQLite database (one-time, ~5-10 min)...\n")
    read.names.sql(names_dmp, sqlFile = tax_sql)
    read.nodes.sql(nodes_dmp, sqlFile = tax_sql)
    cat("  Done:", tax_sql, "\n")
  } else {
    cat("  Using existing taxonomy database:", tax_sql, "\n")
  }

  # Collect unique taxids from BOTH sources together — one getTaxonomy() call
  # serves both, regardless of which taxonomy step a given taxid came from.
  unique_taxids <- unique(na.omit(c(mmseqs_all$lca_taxid, myco_all$myco_lca_taxid)))
  cat("  Unique LCA taxids to expand (UniRef90 + Mycocosm/Phytozome):", length(unique_taxids), "\n")

  tax_lookup <- tryCatch(
    getTaxonomy(as.character(unique_taxids), tax_sql),
    error = function(e) { warning("taxonomizr lookup failed: ", conditionMessage(e)); NULL }
  )

  if (!is.null(tax_lookup)) {
    tax_lookup_dt <- as.data.table(tax_lookup, keep.rownames = "lca_taxid")
    tax_lookup_dt[, lca_taxid := as.integer(lca_taxid)]

    # Rename only ranks present in the returned object; fill missing with NA.
    # getTaxonomy() column names vary across taxonomizr versions.
    present_ranks <- intersect(tax_ranks, names(tax_lookup_dt))
    missing_ranks <- setdiff(tax_ranks, names(tax_lookup_dt))
    if (length(present_ranks) > 0) {
      setnames(tax_lookup_dt, present_ranks, paste0("tax_", present_ranks))
    }
    for (r in missing_ranks) {
      tax_lookup_dt[, (paste0("tax_", r)) := NA_character_]
    }
    if (length(missing_ranks) > 0) {
      cat("  NOTE: ranks not returned by getTaxonomy():", paste(missing_ranks, collapse = ", "), "\n")
    }

    # UniRef90 expansion (existing behavior, tax_* columns)
    lca_table <- mmseqs_all[
      !is.na(lca_taxid)
    ][tax_lookup_dt, on = "lca_taxid", nomatch = NA]

    cat("  UniRef90 rank expansion complete for", nrow(lca_table), "genes.\n")

    if ("tax_superkingdom" %in% names(lca_table)) {
      kingdom_tbl <- lca_table[, .N, by = tax_superkingdom][order(-N)]
      cat("  Kingdom breakdown (UniRef90):\n")
      print(kingdom_tbl)
    }

    # Mycocosm/Phytozome expansion — reuse the same tax_lookup_dt, but
    # re-prefix the rank columns myco_tax_* so they don't collide with the
    # UniRef90 tax_* columns once both are joined into base_dt (Section 5).
    # Join key is renamed to myco_lca_taxid to match myco_all's own columns.
    myco_tax_lookup_dt <- copy(tax_lookup_dt)
    setnames(myco_tax_lookup_dt,
             c("lca_taxid", paste0("tax_", tax_ranks)),
             c("myco_lca_taxid", paste0("myco_tax_", tax_ranks)))

    myco_lca_table <- myco_all[
      !is.na(myco_lca_taxid)
    ][myco_tax_lookup_dt, on = "myco_lca_taxid", nomatch = NA]

    # Informational flags derived from the lineage string — purely additive,
    # no filtering. Subphylum name list is provisional: finalize once the
    # actual lineage strings from the user's genome panel can be inspected,
    # since Mycocosm/NCBI subphylum naming doesn't always align 1:1.
    myco_lca_table[, myco_is_plant := grepl("Viridiplantae", myco_lineage)]
    myco_lca_table[, myco_fungal_subphylum := str_extract(
      myco_lineage,
      "Agaricomycotina|Pezizomycotina|Ustilaginomycotina|Saccharomycotina|Taphrinomycotina"
    )]

    cat("  Mycocosm/Phytozome rank expansion complete for", nrow(myco_lca_table), "genes.\n")
  }
}
cat("\n")

# =============================================================================
# SECTION 5: Build master annotation table
# =============================================================================
# One row per gene. All annotation layers joined by gene_id.
# Genes missing a layer get NAs for that layer's columns.

cat("--- Section 5: Building master annotation table ---\n")

# Base: all genes observed in featureCounts (or union of all annotation files)
all_genes <- unique(c(
  counts_all$gene_id,
  kofam_all$gene_id,
  dbcan_all$gene_id,
  phi_all$gene_id,
  mmseqs_all$gene_id,
  pfam_all$gene_id,
  myco_all$gene_id
))
base_dt <- data.table(gene_id = all_genes)

# Gene-level metadata from featureCounts (contig coords)
gene_meta <- unique(counts_all[, .(gene_id, sample, contig, start, end, strand, length)])
base_dt <- gene_meta[base_dt, on = "gene_id"]

# KOfam: collapse multiple KOs per gene to semicolon-separated string
kofam_collapsed <- kofam_all[
  , .(KO             = paste(unique(KO), collapse = ";"),
      KO_definitions = paste(unique(KO_definition), collapse = ";")
  ), by = gene_id
]
base_dt <- kofam_collapsed[base_dt, on = "gene_id"]

# dbCAN: filter to high-confidence calls (≥ min_tools), collapse families.
# CAZy_substrates is collapsed within this same min_tools-filtered subset, so
# it stays consistent with CAZy_families — a gene with a dbCAN_sub substrate
# hit but family calls below min_tools will show NA here, not a substrate
# with no corresponding family. dbCAN_sub doesn't count toward CAZy_n_tools
# (see 04_dbcan.sh) — it's reported, not voted, regardless of this filter.
dbcan_hicof <- dbcan_all[
  !is.na(CAZy_n_tools) & CAZy_n_tools >= min_tools
][
  , .(CAZy_families  = paste(unique(na.omit(CAZy_HMMER)), collapse = ";"),
      CAZy_substrates = paste(unique(na.omit(CAZy_substrate)), collapse = ";"),
      CAZy_n_tools   = max(CAZy_n_tools, na.rm = TRUE)
  ), by = gene_id
]
# For genes below threshold, keep them with a flag
dbcan_any <- dbcan_all[, .(gene_id, CAZy_HMMER, CAZy_DIAMOND, CAZy_n_tools)]
dbcan_any <- unique(dbcan_any)

base_dt <- dbcan_hicof[base_dt, on = "gene_id"]

# PHI-base: best hit per gene. Filter NA bitscores first so they cannot be
# selected as "best" when order(-phi_bitscore) is applied.
phi_best <- phi_all[!is.na(phi_bitscore)][order(-phi_bitscore)][, .SD[1], by = gene_id]
base_dt <- phi_best[base_dt, on = "gene_id"]

# MMseqs2 LCA taxonomy (gene_id, lca_taxid, lca_rank, lca_name, lineage,
# tax_superkingdom … tax_species). Genes with no hit or unclassified get NAs.
if (!is.null(lca_table)) {
  base_dt <- lca_table[base_dt, on = "gene_id"]
}

# Pfam: collapse multiple domain hits per gene to semicolon-separated strings
# (same convention as KO/CAZy multi-hit genes above).
pfam_collapsed <- pfam_all[
  , .(Pfam_accessions = paste(unique(Pfam_accession), collapse = ";"),
      Pfam_names      = paste(unique(Pfam_name), collapse = ";")
  ), by = gene_id
]
base_dt <- pfam_collapsed[base_dt, on = "gene_id"]

# Mycocosm/Phytozome LCA taxonomy (gene_id, myco_lca_taxid, myco_lca_rank,
# myco_lca_name, myco_lineage, myco_tax_superkingdom … myco_tax_species,
# myco_is_plant, myco_fungal_subphylum). Additive — purely informational,
# no filtering applied. Genes with no hit get NAs for these columns.
if (!is.null(myco_lca_table)) {
  base_dt <- myco_lca_table[base_dt, on = "gene_id"]
}

cat("  Master annotation table:", nrow(base_dt), "genes,",
    ncol(base_dt), "columns\n\n")

# =============================================================================
# SECTION 6: Family-level count matrices
# =============================================================================
# These are the primary inputs for DESeq2 / vegan.
# For each family system (KO, CAZy, PHI phenotype), we:
#   1. Join gene counts with the annotation
#   2. Sum counts per family per sample

cat("--- Section 6: Building family-level count matrices ---\n")

# Helper: build a family × sample count matrix from a long annotation table
# ann_dt must have columns: gene_id, family, sample, count
build_family_matrix <- function(gene_counts, annotation_dt, family_col) {
  # Merge counts with annotation — annotation may have multiple rows per gene
  merged <- annotation_dt[gene_counts[, .(gene_id, sample, count)], on = "gene_id", nomatch = 0]
  merged <- merged[!is.na(get(family_col)) & get(family_col) != ""]

  # Sum counts by family and sample
  summed <- merged[
    , .(count = sum(count, na.rm = TRUE))
    , by = c(family_col, "sample")
  ]

  # Cast to wide matrix
  mat <- dcast(summed, as.formula(paste(family_col, "~ sample")),
               value.var = "count", fill = 0L)
  mat
}

gene_counts_long <- counts_all[, .(gene_id, sample, count)]

# 6a: KO count matrix
# For multi-KO genes, the gene's count is attributed to EACH KO (intentional —
# a gene encoding an enzyme complex contributes to all its KO annotations).
# If this inflation is a concern, divide count by number of KOs per gene.
cat("  Building KO count matrix...\n")
ko_long <- kofam_all[, .(gene_id, family = KO)]
ko_matrix <- build_family_matrix(gene_counts_long, ko_long, "family")
cat("    KO matrix:", nrow(ko_matrix), "KOs ×", ncol(ko_matrix) - 1, "samples\n")

# 6b: CAZy family count matrix (high-confidence calls only)
cat("  Building CAZy family count matrix...\n")
# Expand semicolon-separated families (a gene can have multiple CAZy families).
# gene_id must NOT appear in both .() and by= — the by= clause already adds it
# to the result; including it in .() causes a duplicate column error.
cazy_long <- dbcan_hicof[
  !is.na(CAZy_families) & CAZy_families != "",
  .(family = unlist(strsplit(CAZy_families, ";", fixed = TRUE))),
  by = gene_id
][family != "" & !is.na(family)]
cazy_matrix <- build_family_matrix(gene_counts_long, cazy_long, "family")
cat("    CAZy matrix:", nrow(cazy_matrix), "families ×", ncol(cazy_matrix) - 1, "samples\n")

# 6c: PHI-base phenotype count matrix
cat("  Building PHI-base phenotype count matrix...\n")
phi_long <- phi_all[!is.na(phi_phenotype) & phi_phenotype != "", .(gene_id, family = phi_phenotype)]
phi_matrix <- if (nrow(phi_long) > 0) {
  build_family_matrix(gene_counts_long, phi_long, "family")
} else {
  data.table(family = character())
}
cat("    PHI phenotype matrix:", nrow(phi_matrix), "phenotypes ×",
    max(0, ncol(phi_matrix) - 1), "samples\n\n")

# 6d: Pfam family count matrix
# build_family_matrix() is fully generic over any (gene_id, family) table —
# no new helper logic needed, same instantiation pattern as KO/CAZy/PHI above.
cat("  Building Pfam family count matrix...\n")
pfam_long <- pfam_all[, .(gene_id, family = Pfam_accession)]
pfam_matrix <- if (nrow(pfam_long) > 0) {
  build_family_matrix(gene_counts_long, pfam_long, "family")
} else {
  data.table(family = character())
}
cat("    Pfam matrix:", nrow(pfam_matrix), "families ×",
    max(0, ncol(pfam_matrix) - 1), "samples\n\n")

# =============================================================================
# SECTION 6b: Pseudo-genome (additive log-ratio) normalization
# =============================================================================
# Replicates the normalization approach from the source methodology paper:
# Asparaginase (Pfam PF01112) is treated as a near-single-copy marker, used
# as a proxy for genome equivalents represented by each sample's sequencing
# depth. All four family matrices (KO, CAZy, PHI, Pfam) are normalized against
# it, then log-transformed — mathematically an additive log-ratio (ALR)
# transform using Asparaginase as the reference component.
#
# Why all four matrices, not just Pfam: KO/CAZy/PHI/Pfam matrices are all
# built from the same underlying featureCounts read-pair counts
# (gene_counts_long), just bucketed by different functional-family
# assignments — so numerator and denominator are always on identical units
# regardless of which annotation tool produced the family label. Asparaginase
# acts as a sample-level scaling factor (how much genomic content this
# sample's sequencing represents), independent of which functional system is
# being normalized — the same logic behind single-copy-marker "genome
# equivalents" methods (e.g. MicrobeCensus). NOTE: this is a methodological
# extension beyond the cited paper's literal scope — the paper only had
# Pfam-derived functional counts to normalize, so it never tested this logic
# against KO/CAZy/PHI-style counts from other annotation tools. See README.
#
# Zero-handling: an additive pseudocount is applied to both numerator and
# denominator before the ratio (simpler than Aitchison-style multiplicative
# zero replacement, avoids a new dependency). pseudocount is exposed as a CLI
# flag since the "right" value is somewhat arbitrary and depth-sensitive —
# treat the default as a starting point requiring a sensitivity check, not a
# fixed constant.

cat("--- Section 6b: Pseudo-genome (ALR) normalization ---\n")

# Asparaginase pseudo-genome count per sample: pull the PF01112 row out of
# pfam_matrix. build_family_matrix()'s dcast(..., fill = 0L) already
# guarantees 0 (not NA) for samples with no PF01112 hit.
asparaginase_counts <- setNames(rep(0, length(samples)), samples)
if (nrow(pfam_matrix) > 0 && "PF01112" %in% pfam_matrix$family) {
  pf_row <- pfam_matrix[family == "PF01112"]
  for (s in samples) {
    if (s %in% names(pf_row)) asparaginase_counts[s] <- pf_row[[s]]
  }
} else {
  warning("No PF01112 (Asparaginase) hits found in any sample — pseudo-genome ",
          "normalization will be dominated by the pseudocount for all samples. ",
          "Check that the Pfam step (03b_pfam.sh) ran successfully.")
}

asparaginase_dt <- data.table(sample = samples,
                               asparaginase_pseudogenome_count = asparaginase_counts[samples])
near_zero <- asparaginase_dt[asparaginase_pseudogenome_count <= pseudocount]
if (nrow(near_zero) > 0) {
  cat("  WARNING: samples with near-zero Asparaginase pseudo-genome count",
      "(ALR normalization unreliable for these):\n")
  print(near_zero)
}

# alr_normalize(): for each sample column, ratio = (count + pc) / (denom + pc),
# then natural-log-transform. Applies to any family x sample matrix produced
# by build_family_matrix().
alr_normalize <- function(family_matrix, denom_vec, pseudocount = 1) {
  if (nrow(family_matrix) == 0) return(family_matrix)
  family_col <- names(family_matrix)[1]
  sample_cols <- setdiff(names(family_matrix), family_col)
  result <- copy(family_matrix)
  for (s in sample_cols) {
    denom <- if (s %in% names(denom_vec)) denom_vec[[s]] else NA_real_
    result[[s]] <- log((family_matrix[[s]] + pseudocount) / (denom + pseudocount))
  }
  result
}

denom_vec <- asparaginase_counts

cat("  Normalizing KO, CAZy, PHI, and Pfam matrices against Asparaginase (PF01112)...\n")
ko_matrix_normalized   <- alr_normalize(ko_matrix,   denom_vec, pseudocount)
cazy_matrix_normalized <- alr_normalize(cazy_matrix, denom_vec, pseudocount)
phi_matrix_normalized  <- alr_normalize(phi_matrix,  denom_vec, pseudocount)
pfam_matrix_normalized <- alr_normalize(pfam_matrix, denom_vec, pseudocount)

# Expected transform artifact: normalizing pfam_matrix by its own PF01112 row
# makes that row's ALR value log(1) = 0 for every sample by construction.
# This is NOT a bug and NOT a biological signal ("Asparaginase abundance is
# flat") — see README for this caveat.
cat("\n")

# =============================================================================
# SECTION 6c: OrthoDB genome-equivalent normalization (additive alternative
# to Section 6b's Asparaginase method)
# =============================================================================
# Replicates the normalization approach from "Ectomycorrhizal fungal decay
# traits along a soil nitrogen gradient": instead of one marker gene, reads
# are mapped (01b_orthodb_genecount.sh, pre-assembly) against ~1000-1500
# near-single-copy Dikaryotic orthologs from OrthoDB v12. For each ortholog
# group, read counts are averaged across the group's member sequences (to
# avoid double-counting from within-group sequence redundancy — a read
# homologous to a group can hit multiple per-species member sequences).
# The geometric mean of these per-group averages across all groups is the
# genome-equivalent denominator — far more robust than a single marker gene,
# since any one ortholog's noise/absence is washed out by averaging over
# many. This is ADDITIVE alongside Section 6b's Asparaginase normalization,
# not a replacement — both denominators are kept side by side (see README
# for guidance on which to prefer).
#
# Inputs (paths derived from existing --db-dir/--annotation-dir, no new CLI
# flags needed):
#   {db_dir}/orthodb/gene2og.tsv        sequence_id -> ortholog_group_id,
#                                        built once in 00_setup_databases.sh
#   {ann_dir}/orthodb_genecount/{S}_hits.tsv  per-sample DIAMOND blastx hits
#                                        (R1+R2 merged) vs the OrthoDB DB

cat("--- Section 6c: OrthoDB genome-equivalent normalization ---\n")

gene2og_file <- file.path(db_dir, "orthodb", "gene2og.tsv")
orthodb_dt   <- NULL

if (!file.exists(gene2og_file)) {
  warning("OrthoDB gene2og.tsv not found at ", gene2og_file, " — skipping ",
          "OrthoDB normalization (additive/optional; Asparaginase ",
          "normalization in Section 6b is unaffected).")
} else {
  gene2og <- fread(gene2og_file, header = FALSE, col.names = c("og_seq_id", "og_id"))
  # Group sizes computed once — needed to average each group's total hits
  # across ALL its member sequences, including those with zero hits in a
  # given sample (not just the ones that happened to get a hit).
  og_sizes <- gene2og[, .(n_members = .N), by = og_id]

  geometric_mean_se <- function(x, pseudocount = 1) {
    log_x <- log(x + pseudocount)
    list(geo_mean = exp(mean(log_x)),
         se_log   = sd(log_x) / sqrt(length(log_x)),
         n         = length(log_x))
  }

  orthodb_rows <- lapply(samples, function(s) {
    hits_file <- file.path(ann_dir, "orthodb_genecount", paste0(s, "_hits.tsv"))
    if (!file.exists(hits_file) || file.info(hits_file)$size == 0) {
      return(data.table(sample = s, orthodb_geo_mean = NA_real_,
                        orthodb_se_log = NA_real_, orthodb_n_orthologs = 0L))
    }
    # DIAMOND outfmt 6, no header: qseqid, sseqid, ... (only qseqid/sseqid needed)
    # fill = Inf: never stop early on malformed lines (e.g. a missing newline
    # between R1 and R2 outputs can create a 13-field line); select = c(1,2)
    # means only qseqid/sseqid are loaded regardless of extra fields.
    hits <- fread(hits_file, header = FALSE, sep = "\t", quote = "", fill = Inf,
                  select = c(1, 2), col.names = c("read_id", "og_seq_id"))
    if (nrow(hits) == 0) {
      return(data.table(sample = s, orthodb_geo_mean = NA_real_,
                        orthodb_se_log = NA_real_, orthodb_n_orthologs = 0L))
    }

    # Hits per OrthoDB sequence, then roll up to per-group totals via gene2og
    # (left join so groups with zero hits in this sample still contribute a
    # true zero to the per-group average, not an excluded/missing value).
    hits_per_seq <- hits[, .(n_hits = .N), by = og_seq_id]
    group_totals <- gene2og[hits_per_seq, on = "og_seq_id"][
      , .(total_hits = sum(n_hits, na.rm = TRUE)), by = og_id]
    # Join direction matters: og_sizes must control the row count (the FULL
    # universe of ~1000-1500 OGs), with total_hits looked up and defaulting
    # to NA (then 0 below) for OGs with no hits in this sample. Getting this
    # backwards would silently drop zero-hit OGs from the average.
    group_totals <- group_totals[og_sizes, on = "og_id"]
    group_totals[is.na(total_hits), total_hits := 0]
    group_totals[, group_avg := total_hits / n_members]

    stats <- geometric_mean_se(group_totals$group_avg, pseudocount)
    data.table(sample = s, orthodb_geo_mean = stats$geo_mean,
              orthodb_se_log = stats$se_log, orthodb_n_orthologs = stats$n)
  })

  orthodb_dt <- rbindlist(orthodb_rows)

  near_zero_orthodb <- orthodb_dt[is.na(orthodb_geo_mean) | orthodb_geo_mean <= pseudocount]
  if (nrow(near_zero_orthodb) > 0) {
    cat("  WARNING: samples with missing/near-zero OrthoDB genome-equivalent",
        "estimate (normalization unreliable for these):\n")
    print(near_zero_orthodb)
  }

  orthodb_denom_vec <- setNames(orthodb_dt$orthodb_geo_mean, orthodb_dt$sample)

  cat("  Normalizing KO, CAZy, PHI, and Pfam matrices against OrthoDB genome equivalents...\n")
  ko_matrix_normalized_orthodb   <- alr_normalize(ko_matrix,   orthodb_denom_vec, pseudocount)
  cazy_matrix_normalized_orthodb <- alr_normalize(cazy_matrix, orthodb_denom_vec, pseudocount)
  phi_matrix_normalized_orthodb  <- alr_normalize(phi_matrix,  orthodb_denom_vec, pseudocount)
  pfam_matrix_normalized_orthodb <- alr_normalize(pfam_matrix, orthodb_denom_vec, pseudocount)
}
cat("\n")

# =============================================================================
# SECTION 6d: CAZy direct read mapping — competitive all-kingdom with taxonomy
# Following Bahram 2018 (Nature). Maps raw QC'd reads against CAZy (all
# kingdoms); staxids from DIAMOND hits are expanded to kingdom using
# taxonomizr (already loaded). Produces per-kingdom count matrices and a
# kingdom summary. Normalized by OrthoDB genome equivalents from Section 6c.
# =============================================================================

cat("--- Section 6d: CAZy direct read mapping (Bahram approach) ---\n")

cazy_readmap_dir <- file.path(ann_dir, "cazy_readmap")
cazy_readmap_matrix          <- NULL
cazy_readmap_fungi_matrix    <- NULL
cazy_readmap_bacteria_matrix <- NULL
cazy_readmap_kingdom_dt      <- NULL

if (!dir.exists(cazy_readmap_dir) ||
    length(list.files(cazy_readmap_dir, pattern = "_hits\\.tsv$")) == 0) {
  cat("  [SKIP] No cazy_readmap/ hits files found — run 01c_cazy_readmap.sh first\n")
} else {
  # DIAMOND's staxids column is empty when accession2taxid matching fails for
  # non-standard seqid formats like ACCESSION|FAMILY. We do taxonomy lookup in R:
  #   - Format 1 (GenBank): ACCESSION.ver|FAMILY — extract accession, join taxonmap
  #   - Format 2 (JGI MycoCosm): FAMILY|TAXID|PROTEIN_ID — taxid is field 2
  # Load the cazy_taxonmap.tsv built in 00_setup_databases.sh Section 10.
  taxonmap_file <- file.path(db_dir, "cazy_readmap", "cazy_taxonmap.tsv")
  cazy_taxonmap_dt <- if (file.exists(taxonmap_file)) {
    cat("  Loading CAZy taxonmap for R-side taxonomy lookup...\n")
    fread(taxonmap_file, header = TRUE, sep = "\t", quote = "",
          col.names = c("accession", "taxid"),
          colClasses = c("character", "integer"), key = "accession")
  } else {
    cat("  WARNING: cazy_taxonmap.tsv not found — kingdom assignment will be 'unclassified'\n")
    NULL
  }

  # Reader: load hits, apply Bahram 2018 filters, deduplicate paired reads,
  # extract CAZy family and accession for downstream kingdom assignment.
  read_cazy_readmap <- function(hits_file, pident_min, evalue_max) {
    if (!file.exists(hits_file) || file.info(hits_file)$size == 0) return(NULL)
    # outfmt 6: qseqid(1) sseqid(2) pident(3) length(4) mismatch(5) gapopen(6)
    #           qstart(7) qend(8) sstart(9) send(10) evalue(11) bitscore(12) staxids(13)
    # staxids (col 13) is empty when DIAMOND can't match seqid format — we derive
    # taxids from sseqid directly, so we don't need col 13.
    dt <- fread(hits_file, header = FALSE, sep = "\t", quote = "", fill = Inf,
                select = c(1L, 2L, 3L, 11L, 12L),
                col.names = c("qseqid", "sseqid", "pident", "evalue", "bitscore"))
    if (nrow(dt) == 0) return(NULL)

    # Apply Bahram 2018 final filters (DIAMOND used e-5 at mapping time)
    dt <- dt[pident >= pident_min & evalue <= evalue_max]
    if (nrow(dt) == 0) return(NULL)

    # CAZy database contains three sseqid formats (determined empirically):
    #   GenBank:  ACCESSION.ver|FAMILY[|EC]  e.g. WP_123456.1|GH5  or  AGW.1|GH165|3.2.1.23
    #   JGI simple:   FAMILY|TAXID|PROTEIN   e.g. AA3|452445|Daces1_...
    #   JGI complex:  FAMILY|FAMILY|TAXID|PROTEIN  or  FAMILY|FAMILY|GT5|TAXID|PROTEIN
    # Detect GenBank by presence of '.' in field 1 (versioned accession like WP_123456.1).
    # JGI entries have a CAZy family name as field 1 (no dot).
    # CAZy family pattern: letters + digits + optional underscore-number (subfamily)
    FAMILY_PAT <- "([A-Za-z]+[0-9]+(?:_[0-9]+)?)"
    dt[, field1 := sub("\\|.*$", "", sseqid)]
    dt[, is_genbank := grepl("\\.", field1)]

    dt[, family := {
      f <- rep(NA_character_, .N)
      # GenBank: rightmost field matching family pattern; allow optional trailing
      #   |EC_number (e.g. |3.2.1.23) which trips up a bare $ anchor.
      gb <- is_genbank
      f[gb] <- sub(paste0(".*\\|", FAMILY_PAT, "(?:\\|[^A-Za-z].*)?$"), "\\1", sseqid[gb])
      f[gb & f == sseqid] <- NA_character_
      # JGI: family is always field 1
      jgi <- !is_genbank
      f[jgi] <- sub(paste0("^", FAMILY_PAT, ".*$"), "\\1", sseqid[jgi])
      f[jgi & !grepl("^[A-Za-z]+[0-9]", f)] <- NA_character_
      f
    }]
    dt <- dt[!is.na(family)]
    if (nrow(dt) == 0) return(NULL)

    # Taxid extraction:
    #   GenBank: need taxonmap join — store accession (field1) for joining below
    #   JGI: taxid is the last all-numeric pipe-delimited field in the sseqid
    dt[is_genbank  == TRUE,  accession := field1]
    dt[is_genbank  == FALSE, taxid := {
      parts <- strsplit(sseqid[is_genbank == FALSE], "\\|")
      as.integer(sapply(parts, function(p) {
        num <- p[grepl("^[0-9]+$", p)]
        if (length(num) > 0) tail(num, 1L) else NA_character_
      }))
    }]
    dt[, c("field1", "is_genbank") := NULL]

    # Deduplicate paired reads: R1 and R2 from the same insert can both hit
    # the same family. Strip direction suffixes (/1 /2 .1 .2) to get a shared
    # read-pair ID, then keep only the best-bitscore hit per (pair, family).
    dt[, pair_id := sub("[/. ][12]$", "", qseqid)]
    dt <- dt[dt[, .I[which.max(bitscore)], by = .(pair_id, family)]$V1]

    dt[, .(pair_id, family, taxid, accession, bitscore)]
  }

  # Load hits for all samples that have output files
  cazy_readmap_list <- lapply(samples, function(s) {
    f <- file.path(cazy_readmap_dir, paste0(s, "_hits.tsv"))
    read_cazy_readmap(f, cazy_pident, cazy_evalue)
  })
  names(cazy_readmap_list) <- samples

  # Join GenBank accessions → taxids using cazy_taxonmap_dt (loaded above)
  if (!is.null(cazy_taxonmap_dt)) {
    cazy_readmap_list <- lapply(cazy_readmap_list, function(dt) {
      if (is.null(dt) || nrow(dt) == 0) return(dt)
      needs <- !is.na(dt$accession)
      if (any(needs)) {
        idx <- cazy_taxonmap_dt[dt[needs], on = "accession", which = FALSE]
        dt[needs, taxid := idx$taxid]
      }
      dt[, accession := NULL]
      dt
    })
  } else {
    cazy_readmap_list <- lapply(cazy_readmap_list, function(dt) {
      if (!is.null(dt)) dt[, accession := NULL]
      dt
    })
  }

  n_with_hits <- sum(sapply(cazy_readmap_list, function(x) !is.null(x) && nrow(x) > 0))
  cat("  Samples with CAZy readmap hits:", n_with_hits, "/", length(samples), "\n")

  if (n_with_hits > 0) {
    # Collect all unique taxids across all samples for a single taxonomizr call
    all_taxids <- unique(unlist(lapply(cazy_readmap_list, function(x) {
      if (is.null(x)) return(integer(0))
      x$taxid[!is.na(x$taxid) & x$taxid > 0]
    })))
    cat("  Expanding", length(all_taxids), "unique taxids to kingdom via taxonomizr...\n")

    # getTaxonomy returns a matrix: rows = taxids, cols = superkingdom/phylum/...
    # taxonomizr uses names.dmp + nodes.dmp already loaded from db_dir.
    # The accessionTaxa.sql database path used by getTaxonomy is in db_dir/taxonomy/
    # Use the same taxonomy.sql built by the rest of 08_integrate.R (Section 5).
    # That file is created on first run from names.dmp + nodes.dmp; if it exists
    # we can call getTaxonomy() directly without any additional setup.
    taxa_sql <- file.path(db_dir, "taxonomy", "taxonomy.sql")
    if (file.exists(taxa_sql)) {
      taxon_mat <- getTaxonomy(all_taxids, taxa_sql)
    } else {
      cat("  WARNING: taxonomizr database not found at", taxa_sql, "\n")
      cat("  Run 08_integrate.R once with --annotation-dir and --assembly-dir to build it,\n")
      cat("  or ensure Section 5 (MMseqs2 taxonomy) has run first.\n")
      cat("  Kingdom assignment skipped; all hits will be 'unclassified'\n")
      taxon_mat <- NULL
    }

    # Build taxid → kingdom lookup
    if (!is.null(taxon_mat)) {
      kingdom_lookup <- data.table(
        taxid   = as.integer(rownames(taxon_mat)),
        kingdom = taxon_mat[, "superkingdom"]
      )
      # CAZy fungi superkingdom is "Eukaryota"; distinguish Fungi by checking
      # the "kingdom" rank column where available, else use phylum patterns.
      if ("kingdom" %in% colnames(taxon_mat)) {
        kingdom_lookup[, kingdom := ifelse(!is.na(taxon_mat[, "kingdom"]) &
                                           taxon_mat[, "kingdom"] == "Fungi",
                                           "Fungi", kingdom)]
      }
      kingdom_lookup[is.na(kingdom), kingdom := "unclassified"]
    } else {
      kingdom_lookup <- data.table(taxid = all_taxids, kingdom = "unclassified")
    }

    # Annotate each sample's hits with kingdom, build per-sample family counts
    annotate_kingdom <- function(hits_dt) {
      if (is.null(hits_dt) || nrow(hits_dt) == 0) return(NULL)
      hits_dt[kingdom_lookup, kingdom := i.kingdom, on = "taxid"]
      hits_dt[is.na(kingdom), kingdom := "unclassified"]
      hits_dt
    }

    cazy_readmap_annotated <- lapply(cazy_readmap_list, annotate_kingdom)
    names(cazy_readmap_annotated) <- samples

    # Kingdom summary per sample
    cazy_readmap_kingdom_dt <- rbindlist(lapply(samples, function(s) {
      dt <- cazy_readmap_annotated[[s]]
      if (is.null(dt) || nrow(dt) == 0) {
        return(data.table(sample = s, n_total = 0L, n_fungi = 0L,
                          n_bacteria = 0L, n_unclassified = 0L, pct_fungi = NA_real_))
      }
      n_total  <- nrow(dt)
      n_fungi  <- sum(dt$kingdom == "Fungi", na.rm = TRUE)
      n_bac    <- sum(dt$kingdom == "Bacteria", na.rm = TRUE)
      n_unc    <- sum(dt$kingdom == "unclassified", na.rm = TRUE)
      data.table(sample = s, n_total = n_total, n_fungi = n_fungi,
                 n_bacteria = n_bac, n_unclassified = n_unc,
                 pct_fungi = round(100 * n_fungi / max(n_total, 1), 1))
    }))
    cat("  Kingdom summary (totals across all samples):\n")
    cat("    Total hits :", sum(cazy_readmap_kingdom_dt$n_total), "\n")
    cat("    Fungi      :", sum(cazy_readmap_kingdom_dt$n_fungi), "\n")
    cat("    Bacteria   :", sum(cazy_readmap_kingdom_dt$n_bacteria), "\n")
    cat("    Unclassified:", sum(cazy_readmap_kingdom_dt$n_unclassified), "\n")

    # Build family×sample count matrices per kingdom using existing helper
    make_readmap_input <- function(annotated_list, kingdom_filter) {
      rbindlist(lapply(names(annotated_list), function(s) {
        dt <- annotated_list[[s]]
        if (is.null(dt) || nrow(dt) == 0) return(NULL)
        sub_dt <- if (is.null(kingdom_filter)) dt else dt[kingdom == kingdom_filter]
        if (nrow(sub_dt) == 0) return(NULL)
        sub_dt[, .(sample = s, family, count = 1L)][
          , .(count = sum(count)), by = .(sample, family)]
      }))
    }

    build_readmap_matrix <- function(counts_long) {
      if (is.null(counts_long) || nrow(counts_long) == 0) return(NULL)
      mat <- dcast(counts_long, family ~ sample, value.var = "count", fill = 0L)
      # Add zero columns for samples with no hits
      for (s in samples) {
        if (!s %in% names(mat)) mat[[s]] <- 0L
      }
      setcolorder(mat, c("family", samples))
      mat
    }

    all_long     <- make_readmap_input(cazy_readmap_annotated, NULL)
    fungi_long   <- make_readmap_input(cazy_readmap_annotated, "Fungi")
    bac_long     <- make_readmap_input(cazy_readmap_annotated, "Bacteria")

    cazy_readmap_matrix          <- build_readmap_matrix(all_long)
    cazy_readmap_fungi_matrix    <- build_readmap_matrix(fungi_long)
    cazy_readmap_bacteria_matrix <- build_readmap_matrix(bac_long)

    if (!is.null(cazy_readmap_fungi_matrix))
      cat("  cazy_readmap_fungi_matrix    :", nrow(cazy_readmap_fungi_matrix), "families\n")
    if (!is.null(cazy_readmap_bacteria_matrix))
      cat("  cazy_readmap_bacteria_matrix :", nrow(cazy_readmap_bacteria_matrix), "families\n")
  }
}
cat("\n")

# =============================================================================
# SECTION 7: Per-sample QC summary
# =============================================================================

cat("--- Section 7: Per-sample QC summary ---\n")

summary_rows <- lapply(samples, function(s) {
  n_genes    <- if (!is.null(counts_list[[s]])) nrow(counts_list[[s]])     else 0L
  n_ko       <- if (!is.null(kofam_list[[s]]))  n_distinct(kofam_list[[s]]$gene_id) else 0L
  n_cazy     <- if (!is.null(dbcan_list[[s]])) {
                  sum(dbcan_list[[s]]$CAZy_n_tools >= min_tools, na.rm = TRUE)
                } else 0L
  n_phi      <- if (!is.null(phi_list[[s]]))    nrow(phi_list[[s]])        else 0L
  n_tax <- if (!is.null(mmseqs_list[[s]])) sum(!is.na(mmseqs_list[[s]]$lca_taxid)) else 0L
  n_pfam <- if (!is.null(pfam_list[[s]])) n_distinct(pfam_list[[s]]$gene_id) else 0L
  # NOTE: myco_list elements still use the reader's original column name
  # (lca_taxid) — the myco_lca_taxid rename happens later, only on myco_all.
  n_myco_tax <- if (!is.null(myco_list[[s]])) sum(!is.na(myco_list[[s]]$lca_taxid)) else 0L
  total_reads <- if (!is.null(counts_list[[s]])) sum(counts_list[[s]]$count) else 0L
  asparaginase_count <- asparaginase_counts[[s]]

  data.table(
    sample                  = s,
    genes_predicted         = n_genes,
    total_read_count        = total_reads,
    genes_with_KO           = n_ko,
    pct_KO                  = round(100 * n_ko  / max(n_genes, 1), 1),
    genes_with_CAZy         = n_cazy,
    pct_CAZy                = round(100 * n_cazy / max(n_genes, 1), 1),
    genes_with_PHI_hit      = n_phi,
    pct_PHI                 = round(100 * n_phi  / max(n_genes, 1), 1),
    genes_with_taxonomy     = n_tax,
    pct_taxonomy            = round(100 * n_tax  / max(n_genes, 1), 1),
    genes_with_Pfam         = n_pfam,
    pct_Pfam                = round(100 * n_pfam / max(n_genes, 1), 1),
    genes_with_myco_taxonomy = n_myco_tax,
    pct_myco_taxonomy       = round(100 * n_myco_tax / max(n_genes, 1), 1),
    asparaginase_pseudogenome_count = asparaginase_count
  )
})

summary_dt <- rbindlist(summary_rows)

# Join stand_type from stakes.csv (--metadata) if provided.
# stake_id is extracted from the sample name using the same regex as the QC
# report: GEO-{stake_id}[_Cleaned]_S{N} → stake_id = the NN-NN part.
if (!is.null(metadata_file) && file.exists(metadata_file)) {
  meta <- fread(metadata_file, colClasses = "character")
  # stakes.csv uses column "stake" (not "stake_id"); stand_type is the group label.
  # Regex extracts the NN-NN stake number from names like:
  #   GEO-00-06_S1, GEO-03-13_Cleaned_S8, GEO-19-17-Cleaned
  # The .*$ at the end handles any suffix (_Cleaned, _SN, -Cleaned, etc.)
  if (all(c("stake", "stand_type") %in% names(meta))) {
    summary_dt[, stake_id := sub("^GEO-([0-9]+-[0-9]+).*$", "\\1", sample)]
    summary_dt <- merge(summary_dt, meta[, .(stake_id = stake, stand_type)],
                        by = "stake_id", all.x = TRUE, sort = FALSE)
    setcolorder(summary_dt, c("sample", "stake_id", "stand_type"))
    cat("  Merged stand_type for", sum(!is.na(summary_dt$stand_type)), "of",
        nrow(summary_dt), "samples\n")
  } else {
    cat("  WARNING: --metadata file missing 'stake' or 'stand_type' column — skipping metadata join\n")
    cat("  Columns found:", paste(names(meta), collapse = ", "), "\n")
  }
} else if (!is.null(metadata_file)) {
  cat("  WARNING: --metadata file not found:", metadata_file, "— skipping metadata join\n")
}

# Join OrthoDB genome-equivalent estimates into summary_dt.
if (!is.null(orthodb_dt)) {
  summary_dt <- merge(summary_dt, orthodb_dt, by = "sample", all.x = TRUE, sort = FALSE)
  cat("  Merged OrthoDB genome equivalents for",
      sum(!is.na(summary_dt$orthodb_geo_mean)), "of", nrow(summary_dt), "samples\n")
}

# Join CAZy readmap kingdom counts into summary_dt.
if (!is.null(cazy_readmap_kingdom_dt)) {
  summary_dt <- merge(summary_dt,
                      cazy_readmap_kingdom_dt[, .(sample, n_cazy_readmap_total = n_total,
                                                   n_cazy_readmap_fungi = n_fungi,
                                                   n_cazy_readmap_bacteria = n_bacteria,
                                                   pct_cazy_readmap_fungi = pct_fungi)],
                      by = "sample", all.x = TRUE, sort = FALSE)
  cat("  Merged CAZy readmap kingdom counts for",
      sum(!is.na(summary_dt$n_cazy_readmap_total)), "of", nrow(summary_dt), "samples\n")
}

print(summary_dt)
cat("\n")

# =============================================================================
# SECTION 8: Write outputs
# =============================================================================

cat("--- Section 8: Writing outputs ---\n")

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", na = "NA", quote = FALSE)
  cat("  Wrote:", path, paste0("(", nrow(dt), " rows)\n"))
}

write_tsv(base_dt,       file.path(out_dir, "gene_annotations.tsv"))
write_tsv(count_matrix,  file.path(out_dir, "gene_counts_raw.tsv"))
write_tsv(ko_matrix,     file.path(out_dir, "ko_count_matrix.tsv"))
write_tsv(cazy_matrix,   file.path(out_dir, "cazy_count_matrix.tsv"))
write_tsv(phi_matrix,    file.path(out_dir, "phi_count_matrix.tsv"))
write_tsv(pfam_matrix,   file.path(out_dir, "pfam_count_matrix.tsv"))
write_tsv(summary_dt,    file.path(out_dir, "summary_stats.tsv"))

# Pseudo-genome (ALR) normalized matrices — see Section 6b for the transform
# and the README for the all-four-matrices normalization rationale.
write_tsv(ko_matrix_normalized,   file.path(out_dir, "ko_count_matrix_normalized.tsv"))
write_tsv(cazy_matrix_normalized, file.path(out_dir, "cazy_count_matrix_normalized.tsv"))
write_tsv(phi_matrix_normalized,  file.path(out_dir, "phi_count_matrix_normalized.tsv"))
write_tsv(pfam_matrix_normalized, file.path(out_dir, "pfam_count_matrix_normalized.tsv"))
write_tsv(asparaginase_dt,        file.path(out_dir, "asparaginase_pseudogenome_counts.tsv"))

# OrthoDB-normalized matrices — additive alongside the Asparaginase-normalized
# ones above (Section 6c). Only written if gene2og.tsv / hit tables were
# available; otherwise this is silently skipped (optional/additive method).
if (!is.null(orthodb_dt)) {
  write_tsv(ko_matrix_normalized_orthodb,   file.path(out_dir, "ko_count_matrix_normalized_orthodb.tsv"))
  write_tsv(cazy_matrix_normalized_orthodb, file.path(out_dir, "cazy_count_matrix_normalized_orthodb.tsv"))
  write_tsv(phi_matrix_normalized_orthodb,  file.path(out_dir, "phi_count_matrix_normalized_orthodb.tsv"))
  write_tsv(pfam_matrix_normalized_orthodb, file.path(out_dir, "pfam_count_matrix_normalized_orthodb.tsv"))
  write_tsv(orthodb_dt,                     file.path(out_dir, "orthodb_genome_equivalents.tsv"))
}

# CAZy readmap matrices (Section 6d) — competitive all-kingdom read mapping
# following Bahram 2018. Only written if 01c_cazy_readmap.sh has been run.
if (!is.null(cazy_readmap_kingdom_dt)) {
  write_tsv(cazy_readmap_kingdom_dt, file.path(out_dir, "cazy_readmap_kingdom_summary.tsv"))
}
if (!is.null(cazy_readmap_matrix)) {
  write_tsv(cazy_readmap_matrix, file.path(out_dir, "cazy_readmap_matrix.tsv"))
}
if (!is.null(cazy_readmap_fungi_matrix)) {
  write_tsv(cazy_readmap_fungi_matrix, file.path(out_dir, "cazy_readmap_fungi_matrix.tsv"))
  if (!is.null(orthodb_dt)) {
    cazy_readmap_fungi_normalized <- alr_normalize(cazy_readmap_fungi_matrix,
                                                   orthodb_denom_vec, pseudocount)
    write_tsv(cazy_readmap_fungi_normalized,
              file.path(out_dir, "cazy_readmap_fungi_matrix_normalized_orthodb.tsv"))
  }
}
if (!is.null(cazy_readmap_bacteria_matrix)) {
  write_tsv(cazy_readmap_bacteria_matrix, file.path(out_dir, "cazy_readmap_bacteria_matrix.tsv"))
}

# Also save as R objects for direct use in downstream analysis
save_objects <- c("base_dt", "count_matrix", "ko_matrix", "cazy_matrix", "phi_matrix",
                  "pfam_matrix", "ko_matrix_normalized", "cazy_matrix_normalized",
                  "phi_matrix_normalized", "pfam_matrix_normalized",
                  "asparaginase_dt", "summary_dt", "orthodb_dt")
if (!is.null(orthodb_dt)) {
  save_objects <- c(save_objects, "ko_matrix_normalized_orthodb", "cazy_matrix_normalized_orthodb",
                    "phi_matrix_normalized_orthodb", "pfam_matrix_normalized_orthodb")
}
if (!is.null(cazy_readmap_fungi_matrix)) {
  save_objects <- c(save_objects, "cazy_readmap_matrix", "cazy_readmap_fungi_matrix",
                    "cazy_readmap_bacteria_matrix", "cazy_readmap_kingdom_dt")
}
save(list = save_objects, file = file.path(out_dir, "integrated_data.RData"))
cat("  Wrote:", file.path(out_dir, "integrated_data.RData"), "(all objects)\n")

cat("\n============================================================\n")
cat("Integration complete:", format(Sys.time()), "\n")
cat("Outputs in:", out_dir, "\n")
cat("\nTo load in R:\n")
cat("  load('", file.path(out_dir, "integrated_data.RData"), "')\n", sep = "")
cat("\nKey analysis tables:\n")
cat("  ko_matrix / ko_matrix_normalized       — KO × sample (raw / ALR vs Asparaginase) → DESeq2 / vegan\n")
cat("  cazy_matrix / cazy_matrix_normalized   — CAZy family × sample (raw / ALR vs Asparaginase)\n")
cat("  phi_matrix / phi_matrix_normalized     — PHI phenotype × sample (raw / ALR vs Asparaginase)\n")
cat("  pfam_matrix / pfam_matrix_normalized   — Pfam family × sample (raw / ALR vs Asparaginase)\n")
cat("  *_matrix_normalized_orthodb            — same four matrices, ALR vs OrthoDB genome equivalents (if available)\n")
cat("  asparaginase_dt                        — per-sample ALR denominator (PF01112 reads)\n")
cat("  orthodb_dt                             — per-sample ALR denominator (OrthoDB geometric mean + SE, if available)\n")
cat("  base_dt                                — full gene annotation table  → custom queries\n")
cat("  NOTE: *_normalized matrices are already log-transformed — do not feed\n")
cat("        them into DESeq2/vegan, which expect raw counts.\n")
cat("============================================================\n")
