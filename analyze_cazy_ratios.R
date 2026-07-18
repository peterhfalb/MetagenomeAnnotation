#!/usr/bin/env Rscript
# analyze_cazy_ratios.R
# Compute fungal:bacterial CAZyme log-ratios from cazy_readmap matrices at four
# levels: overall (per sample), per family, per substrate, and per origin.
#
# Inputs:
#   cazy_readmap_fungi_matrix.tsv   }
#   cazy_readmap_bacteria_matrix.tsv} from 08_integrate.R --out-dir
#   Cazy_lookup_table_fungi.csv       maps families → substrate + origin
#   stakes.csv                        metadata: stake (sample ID), stand_type
#
# Usage:
#   Rscript analyze_cazy_ratios.R \
#     --integrated-dir /path/to/integrated \
#     --lookup-table   /path/to/Cazy_lookup_table_fungi.csv \
#     --metadata       /path/to/stakes.csv \
#     [--out-dir       /path/to/output] \
#     [--pseudocount   1.0]

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

option_list <- list(
  make_option("--integrated-dir", type = "character", default = NULL,
              help = "Path to 08_integrate.R output directory containing readmap matrices (required)"),
  make_option("--lookup-table",   type = "character", default = NULL,
              help = "Path to Cazy_lookup_table_fungi.csv (required)"),
  make_option("--metadata",       type = "character", default = NULL,
              help = "Path to stakes.csv with columns: stake, stand_type (required)"),
  make_option("--out-dir",        type = "character", default = NULL,
              help = "Output directory [default: {integrated-dir}/cazy_ratios/]"),
  make_option("--pseudocount",    type = "double",    default = 1.0,
              help = "Pseudocount added before log-ratio [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt[["integrated-dir"]]) || is.null(opt[["lookup-table"]]) ||
    is.null(opt[["metadata"]])) {
  stop("--integrated-dir, --lookup-table, and --metadata are required.", call. = FALSE)
}

integrated_dir <- opt[["integrated-dir"]]
lookup_file    <- opt[["lookup-table"]]
metadata_file  <- opt[["metadata"]]
out_dir        <- if (is.null(opt[["out-dir"]])) file.path(integrated_dir, "cazy_ratios") else opt[["out-dir"]]
pseudocount    <- opt[["pseudocount"]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("CAZy fungi:bacteria ratio analysis\n")
cat("Integrated dir :", integrated_dir, "\n")
cat("Lookup table   :", lookup_file, "\n")
cat("Output dir     :", out_dir, "\n")
cat("Pseudocount    :", pseudocount, "\n")
cat("============================================================\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", quote = FALSE, na = "NA")
  cat("  Wrote:", path, "\n")
}

# Wilcoxon summary across stand_type groups for a given grouping column.
# Returns one row per group with n, median, and p-value for each stand_type pair.
summarize_ratio <- function(dt, group_col) {
  stand_types <- sort(unique(dt$stand_type))
  groups      <- sort(unique(dt[[group_col]]))

  rbindlist(lapply(groups, function(g) {
    sub <- dt[get(group_col) == g]
    grp_stats <- lapply(stand_types, function(st) {
      v <- sub[stand_type == st, log_ratio]
      list(n = length(v), median = median(v, na.rm = TRUE))
    })
    names(grp_stats) <- stand_types

    # Wilcoxon between first two stand_types (AM vs. ECM or whatever pair exists)
    w <- NULL
    if (length(stand_types) == 2) {
      v1 <- sub[stand_type == stand_types[1], log_ratio]
      v2 <- sub[stand_type == stand_types[2], log_ratio]
      if (length(v1) >= 2 && length(v2) >= 2) {
        w <- tryCatch(wilcox.test(v1, v2), error = function(e) NULL)
      }
    }

    row <- data.table(group = g)
    for (st in stand_types) {
      row[[paste0("n_",      st)]] <- grp_stats[[st]]$n
      row[[paste0("median_", st)]] <- grp_stats[[st]]$median
    }
    if (length(stand_types) == 2) {
      row[, delta_median := grp_stats[[stand_types[2]]]$median - grp_stats[[stand_types[1]]]$median]
    }
    row[, wilcox_p := if (!is.null(w)) w$p.value else NA_real_]
    row
  }))
}

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("--- Loading matrices ---\n")
fung_file <- file.path(integrated_dir, "cazy_readmap_fungi_matrix.tsv")
bact_file <- file.path(integrated_dir, "cazy_readmap_bacteria_matrix.tsv")

if (!file.exists(fung_file)) stop("Fungi matrix not found: ", fung_file, call. = FALSE)
if (!file.exists(bact_file)) stop("Bacteria matrix not found: ", bact_file, call. = FALSE)

fung <- fread(fung_file)
bact <- fread(bact_file)
cat("  Fungi matrix  :", nrow(fung), "families ×", ncol(fung) - 1, "samples\n")
cat("  Bact  matrix  :", nrow(bact), "families ×", ncol(bact) - 1, "samples\n")

# Melt to long format
fung_long <- melt(fung, id.vars = "family", variable.name = "sample",
                  value.name = "n_fungi",    variable.factor = FALSE)
bact_long <- melt(bact, id.vars = "family", variable.name = "sample",
                  value.name = "n_bacteria", variable.factor = FALSE)

# Full outer join — families absent from one matrix get 0, not NA
hits <- merge(fung_long, bact_long, by = c("family", "sample"), all = TRUE)
hits[is.na(n_fungi),    n_fungi    := 0L]
hits[is.na(n_bacteria), n_bacteria := 0L]
hits[, log_ratio := log((n_fungi + pseudocount) / (n_bacteria + pseudocount))]

cat("  Combined hits :", nrow(hits), "family × sample rows\n\n")

# Load metadata
cat("--- Loading metadata ---\n")
meta <- fread(metadata_file)
if (!"stake"      %in% names(meta)) stop("metadata must have a 'stake' column",      call. = FALSE)
if (!"stand_type" %in% names(meta)) stop("metadata must have a 'stand_type' column", call. = FALSE)
meta <- meta[, .(sample = as.character(stake), stand_type = as.character(stand_type))]

n_before <- nrow(hits)
hits <- hits[meta, on = "sample", nomatch = 0]
n_after  <- nrow(hits)
if (n_before != n_after)
  cat("  WARNING:", (n_before - n_after) / length(unique(fung_long$family)),
      "samples dropped (no metadata match)\n")
cat("  Samples matched:", uniqueN(hits$sample), "\n")
cat("  Stand types    :", paste(sort(unique(hits$stand_type)), collapse = ", "), "\n\n")

# Load lookup table
cat("--- Loading lookup table ---\n")
lookup_raw <- fread(lookup_file)
lookup <- lookup_raw[, .(
  family       = as.character(Family),
  substrate    = as.character(Substrate),
  origin       = as.character(Origin),
  decay_signal = as.character(Decay_type_signal)
)]
cat("  Lookup entries :", nrow(lookup), "rows (", uniqueN(lookup$family), "unique families)\n")
cat("  Substrates     :", paste(sort(unique(lookup$substrate)), collapse = ", "), "\n")
cat("  Origins        :", paste(sort(unique(lookup$origin)),    collapse = ", "), "\n\n")

# =============================================================================
# 2. LEVEL 1 — OVERALL PER SAMPLE
# =============================================================================

cat("--- Level 1: overall per sample ---\n")
overall <- hits[, .(
  n_fungi    = sum(n_fungi),
  n_bacteria = sum(n_bacteria),
  log_ratio  = log((sum(n_fungi) + pseudocount) / (sum(n_bacteria) + pseudocount))
), by = .(sample, stand_type)]
setorder(overall, stand_type, sample)

cat("  Rows:", nrow(overall), "\n")
cat("  Median log-ratio by stand_type:\n")
print(overall[, .(median_log_ratio = round(median(log_ratio), 3),
                   n_samples = .N), by = stand_type])

summary_overall <- overall[, .(
  n         = .N,
  median    = median(log_ratio),
  mean      = mean(log_ratio),
  sd        = sd(log_ratio)
), by = stand_type]
cat("\n")

# =============================================================================
# 3. LEVEL 2 — PER FAMILY
# =============================================================================

cat("--- Level 2: per family ---\n")

# Left join lookup onto hits: families not in lookup get NA substrate/origin
by_family <- copy(hits[, .(sample, stand_type, family, n_fungi, n_bacteria, log_ratio)])

# A family can have multiple rows in lookup (multiple substrates) — we want only
# unique substrate/origin per family for the annotation columns here; use the
# first-occurring row per family (deterministic given stable file order).
lookup_primary <- lookup[!duplicated(family)]
by_family <- lookup_primary[by_family, on = "family"]   # left join, NAs for unmapped
setorder(by_family, stand_type, sample, family)

n_mapped <- sum(!is.na(by_family$substrate))
cat("  Family × sample rows :", nrow(by_family), "\n")
cat("  Rows with substrate annotation:", n_mapped, "/", nrow(by_family),
    sprintf("(%.0f%% of families by hits)\n", 100 * n_mapped / nrow(by_family)))
cat("\n")

# =============================================================================
# 4. LEVEL 3 — PER SUBSTRATE
# =============================================================================

cat("--- Level 3: per substrate ---\n")

# Inner join: families not in lookup are excluded from substrate analysis.
# allow.cartesian = TRUE: families with multiple substrate rows are duplicated
# (each copy contributes to its respective substrate group).
hits_sub <- lookup[hits, on = "family", allow.cartesian = TRUE, nomatch = 0]
by_substrate <- hits_sub[, .(
  n_fungi    = sum(n_fungi),
  n_bacteria = sum(n_bacteria),
  n_families = uniqueN(family),
  log_ratio  = log((sum(n_fungi) + pseudocount) / (sum(n_bacteria) + pseudocount))
), by = .(sample, stand_type, substrate)]
setorder(by_substrate, substrate, stand_type, sample)

cat("  Rows:", nrow(by_substrate), "\n")
cat("  Substrates covered:", uniqueN(by_substrate$substrate), "\n")

summary_substrate <- summarize_ratio(by_substrate, "substrate")
setorder(summary_substrate, group)
cat("  Summary:\n")
print(summary_substrate)
cat("\n")

# =============================================================================
# 5. LEVEL 4 — PER ORIGIN
# =============================================================================

cat("--- Level 4: per origin ---\n")

hits_ori <- lookup[hits, on = "family", allow.cartesian = TRUE, nomatch = 0]
by_origin <- hits_ori[, .(
  n_fungi    = sum(n_fungi),
  n_bacteria = sum(n_bacteria),
  n_families = uniqueN(family),
  log_ratio  = log((sum(n_fungi) + pseudocount) / (sum(n_bacteria) + pseudocount))
), by = .(sample, stand_type, origin)]
setorder(by_origin, origin, stand_type, sample)

cat("  Rows:", nrow(by_origin), "\n")

summary_origin <- summarize_ratio(by_origin, "origin")
setorder(summary_origin, group)
cat("  Summary:\n")
print(summary_origin)
cat("\n")

# =============================================================================
# 6. WRITE OUTPUTS
# =============================================================================

cat("--- Writing outputs ---\n")
write_tsv(overall,           file.path(out_dir, "ratio_overall.tsv"))
write_tsv(summary_overall,   file.path(out_dir, "ratio_summary_overall.tsv"))
write_tsv(by_family,         file.path(out_dir, "ratio_by_family.tsv"))
write_tsv(by_substrate,      file.path(out_dir, "ratio_by_substrate.tsv"))
write_tsv(summary_substrate, file.path(out_dir, "ratio_summary_by_substrate.tsv"))
write_tsv(by_origin,         file.path(out_dir, "ratio_by_origin.tsv"))
write_tsv(summary_origin,    file.path(out_dir, "ratio_summary_by_origin.tsv"))

save(overall, summary_overall,
     by_family,
     by_substrate, summary_substrate,
     by_origin,    summary_origin,
     file = file.path(out_dir, "cazy_ratios.RData"))
cat("  Wrote: cazy_ratios.RData\n\n")

cat("============================================================\n")
cat("Done.\n")
cat("  Log-ratio > 0  = fungi-dominated\n")
cat("  Log-ratio < 0  = bacteria-dominated\n")
cat("  Log-ratio = 0  = equal contributions\n")
cat("  Pseudocount    :", pseudocount, "(applied before log)\n")
cat("============================================================\n")
