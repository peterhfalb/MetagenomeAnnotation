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
#
# Interactive use (in an R session):
#   Set variables in the CONFIG block below, then source() this file.
#
# Outputs (all in out-dir/):
#   gene_annotations.tsv    one row per gene — contig coords + all annotation layers
#   gene_counts_raw.tsv     raw count matrix: genes × samples (sparse — one assembly
#                           per sample means each gene appears in one sample only)
#   ko_count_matrix.tsv     KO × sample count matrix       ← main input for DESeq2
#   cazy_count_matrix.tsv   CAZy family × sample matrix    ← main input for DESeq2
#   phi_count_matrix.tsv    PHI-base phenotype × sample    ← pathogenicity analysis
#   summary_stats.tsv       per-sample QC: gene counts, annotation rates
#
# NOTE on count normalization:
#   Outputs are RAW read counts. Normalize in R before differential analysis:
#   - DESeq2: pass raw counts directly (DESeq2 handles its own normalization)
#   - vegan:  normalize with decostand() (e.g., "total", "log", "hellinger")
#   - Do NOT pre-normalize before DESeq2.
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
              help = "LCA: retain hits within this fraction of top bitscore [default: %default]")
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
    gene_id = sub("_mRNA$", "", Geneid),
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
  if (!file.exists(file) || file.info(file)$size == 0) return(data.table(gene_id = character(), KO = character(), KO_definition = character()))
  dt <- fread(file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
  # Columns: threshold_flag, gene_id, KO, threshold, score, evalue, definition
  sig <- dt[V1 == "*"]
  if (nrow(sig) == 0) return(data.table(gene_id = character(), KO = character(), KO_definition = character()))
  sig[, .(gene_id = normalize_gene_id(V2), KO = V3, KO_definition = if (ncol(dt) >= 7) V7 else NA_character_)]
}

# Read dbCAN3 overview.txt for one sample.
# Returns data.table: gene_id, CAZy_HMMER, CAZy_DIAMOND, n_tools
# Strips position suffixes from HMMER calls (e.g. "GH3(1-200)" → "GH3").
read_dbcan <- function(file) {
  if (is.null(file) || !file.exists(file) || file.info(file)$size == 0) return(data.table(gene_id = character()))
  dt <- fread(file, header = TRUE, sep = "\t", quote = "")
  if (nrow(dt) == 0) return(data.table(gene_id = character()))
  # Column names vary slightly by version; normalise
  setnames(dt, c("Gene.ID", "Gene ID"), c("gene_id", "gene_id"), skip_absent = TRUE)
  if (!"gene_id" %in% names(dt)) setnames(dt, 1, "gene_id")
  dt[, gene_id := normalize_gene_id(gene_id)]

  # Identify columns by name patterns — use ignore.case instead of (?i) inline
  # flag so it works with R's default regex engine without requiring perl=TRUE.
  # v3 used "HMMER"; v4 uses "dbCAN_hmm" — both contain "hmm".
  hmmer_col   <- grep("hmm",     names(dt), value = TRUE, ignore.case = TRUE)[1]
  diamond_col <- grep("diamond", names(dt), value = TRUE, ignore.case = TRUE)[1]
  tools_col   <- grep("#.*tool|ofTool|n_tool", names(dt), value = TRUE, ignore.case = TRUE)[1]

  result <- data.table(
    gene_id       = dt[["gene_id"]],
    CAZy_HMMER    = if (!is.na(hmmer_col))   str_remove(dt[[hmmer_col]],   "\\(.*\\)") else NA_character_,
    CAZy_DIAMOND  = if (!is.na(diamond_col)) dt[[diamond_col]] else NA_character_,
    CAZy_n_tools  = if (!is.na(tools_col))   as.integer(dt[[tools_col]])   else NA_integer_
  )
  # Normalise no-hit strings to NA — v3 used "N/A", v4 uses "-"
  result[CAZy_HMMER   %in% c("N/A", "-"), CAZy_HMMER   := NA_character_]
  result[CAZy_DIAMOND %in% c("N/A", "-"), CAZy_DIAMOND := NA_character_]
  result
}

# Parse PHI-base stitle into components.
# PHI-base headers follow: PHI:ACCESSION|gene|pathogen|host|phenotype|...
# or older format: gene_name pathogen_species | phenotype
# Returns a list with: phi_accession, phi_gene, phi_pathogen, phi_phenotype
parse_phi_stitle <- function(stitle) {
  # Try pipe-delimited format (PHI:ID|gene|pathogen|host|phenotype...)
  if (grepl("\\|", stitle)) {
    parts <- str_split(stitle, "\\|")[[1]]
    list(
      phi_accession = str_extract(parts[1], "PHI:[0-9]+"),
      phi_gene      = if (length(parts) >= 2) trimws(parts[2]) else NA_character_,
      phi_pathogen  = if (length(parts) >= 3) trimws(parts[3]) else NA_character_,
      phi_phenotype = if (length(parts) >= 5) trimws(parts[5]) else NA_character_
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

# Normalize MetaEuk gene IDs to the first 4 pipe-delimited fields
# (targetID|contigID|strand|startPos). MetaEuk FASTA headers carry additional
# fields (evalue, nExons, coordinates) that featureCounts never sees — it
# reads gene IDs from the GFF3 Parent= attribute, which only has 4 fields.
# Applying this to all annotation layers ensures the join key matches.
normalize_gene_id <- function(ids) {
  sub("^([^|]+\\|[^|]+\\|[^|]+\\|[^|]+)(\\|.*)?$", "\\1", ids)
}

# Find dbCAN overview file — filename changed across dbCAN v4 sub-versions.
# Searches the sample's output directory for any known candidate name.
find_dbcan_overview <- function(dbcan_dir) {
  candidates <- file.path(dbcan_dir,
    c("overview.txt", "overview.tsv", "CAZyme_annotation.tsv"))
  for (f in candidates) {
    if (file.exists(f) && file.info(f)$size > 0) return(f)
  }
  NULL
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
mmseqs_list  <- list()  # MMseqs2 taxonomy per sample

for (s in samples) {
  cat("  Reading:", s, "\n")

  # featureCounts
  counts_list[[s]] <- read_featurecounts(
    file.path(ann_dir, "featurecounts", paste0(s, "_counts.txt")), s)

  # KOfamScan (reads the full detail-tsv; mapper is derived from it in read_kofam)
  kofam_list[[s]] <- read_kofam(file.path(ann_dir, "kofam", paste0(s, "_kofam.tsv")))

  # dbCAN3 — find actual overview file (name varies across v4 sub-versions)
  dbcan_list[[s]] <- read_dbcan(find_dbcan_overview(file.path(ann_dir, "dbcan", s)))

  # PHI-base
  phi_list[[s]] <- read_phibase(file.path(ann_dir, "phibase", paste0(s, "_phibase.tsv")))

  # MMseqs2 taxonomy (_lca.tsv — LCA already computed by MMseqs2)
  mmseqs_list[[s]] <- read_mmseqs_taxonomy(
    file.path(ann_dir, "mmseqs_taxonomy", paste0(s, "_lca.tsv")))
}

# Combine across samples
counts_all <- rbindlist(counts_list,  use.names = TRUE, fill = TRUE)
kofam_all  <- rbindlist(kofam_list,   use.names = TRUE, fill = TRUE)
dbcan_all  <- rbindlist(dbcan_list,   use.names = TRUE, fill = TRUE)
phi_all    <- rbindlist(phi_list,     use.names = TRUE, fill = TRUE)
mmseqs_all <- rbindlist(mmseqs_list,  use.names = TRUE, fill = TRUE)

cat("\n  Genes with featureCounts data :", nrow(counts_all), "\n")
cat("  KOfam significant hits        :", nrow(kofam_all), "\n")
cat("  dbCAN gene calls              :", nrow(dbcan_all), "\n")
cat("  PHI-base hits                 :", nrow(phi_all), "\n")
cat("  MMseqs2 taxonomy assignments  :", sum(!is.na(mmseqs_all$lca_taxid)), "\n\n")

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

cat("--- Section 4: Expanding MMseqs2 LCA taxon IDs to standard ranks ---\n")

tax_sql   <- file.path(db_dir, "taxonomy", "taxonomy.sql")
names_dmp <- file.path(db_dir, "taxonomy", "names.dmp")
nodes_dmp <- file.path(db_dir, "taxonomy", "nodes.dmp")
tax_ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

lca_table <- NULL

if (!file.exists(names_dmp) || !file.exists(nodes_dmp)) {
  warning("NCBI taxonomy files not found — skipping rank expansion.\n",
          "  Expected: ", names_dmp)
} else if (nrow(mmseqs_all) == 0 || all(is.na(mmseqs_all$lca_taxid))) {
  warning("No classified MMseqs2 hits — skipping rank expansion.")
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

  unique_taxids <- unique(na.omit(mmseqs_all$lca_taxid))
  cat("  Unique LCA taxids to expand:", length(unique_taxids), "\n")

  tax_lookup <- tryCatch(
    getTaxonomy(as.character(unique_taxids), tax_sql),
    error = function(e) { warning("taxonomizr lookup failed: ", conditionMessage(e)); NULL }
  )

  if (!is.null(tax_lookup)) {
    tax_lookup_dt <- as.data.table(tax_lookup, keep.rownames = "lca_taxid")
    tax_lookup_dt[, lca_taxid := as.integer(lca_taxid)]
    setnames(tax_lookup_dt, tax_ranks, paste0("tax_", tax_ranks))

    lca_table <- mmseqs_all[
      !is.na(lca_taxid)
    ][tax_lookup_dt, on = "lca_taxid", nomatch = NA]

    cat("  Rank expansion complete for", nrow(lca_table), "genes.\n")

    kingdom_tbl <- lca_table[, .N, by = tax_superkingdom][order(-N)]
    cat("  Kingdom breakdown:\n")
    print(kingdom_tbl)
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
  mmseqs_all$gene_id
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

# dbCAN: filter to high-confidence calls (≥ min_tools), collapse families
dbcan_hicof <- dbcan_all[
  !is.na(CAZy_n_tools) & CAZy_n_tools >= min_tools
][
  , .(CAZy_families = paste(unique(na.omit(CAZy_HMMER)), collapse = ";"),
      CAZy_n_tools  = max(CAZy_n_tools, na.rm = TRUE)
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
  total_reads <- if (!is.null(counts_list[[s]])) sum(counts_list[[s]]$count) else 0L

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
    pct_taxonomy            = round(100 * n_tax  / max(n_genes, 1), 1)
  )
})

summary_dt <- rbindlist(summary_rows)
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
write_tsv(summary_dt,    file.path(out_dir, "summary_stats.tsv"))

# Also save as R objects for direct use in downstream analysis
save(base_dt, count_matrix, ko_matrix, cazy_matrix, phi_matrix, summary_dt,
     file = file.path(out_dir, "integrated_data.RData"))
cat("  Wrote:", file.path(out_dir, "integrated_data.RData"), "(all objects)\n")

cat("\n============================================================\n")
cat("Integration complete:", format(Sys.time()), "\n")
cat("Outputs in:", out_dir, "\n")
cat("\nTo load in R:\n")
cat("  load('", file.path(out_dir, "integrated_data.RData"), "')\n", sep = "")
cat("\nKey analysis tables:\n")
cat("  ko_matrix    — KO × sample count matrix      → DESeq2 / vegan\n")
cat("  cazy_matrix  — CAZy family × sample matrix   → DESeq2 / vegan\n")
cat("  phi_matrix   — PHI phenotype × sample matrix → DESeq2 / vegan\n")
cat("  base_dt      — full gene annotation table    → custom queries\n")
cat("============================================================\n")
