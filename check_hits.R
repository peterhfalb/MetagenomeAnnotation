#!/usr/bin/env Rscript
# check_hits.R
# Diagnostics on raw CAZy readmap hits (01c_cazy_readmap.sh output). Read-only;
# applies the same pident/e-value filters as 08_integrate.R, then reports:
#   1. Multi-family headers, where the label kept depends on label order
#   2. Unassigned "family 0" bins (GH0, AA0, ...)
#   3. Headers whose first field is neither a CAZy label nor a versioned accession
#   4. Paired-read double counts removed by merging to family before dedup
#
# Usage (on a compute node; hits files can be large):
#   module load R
#   Rscript check_hits.R /path/to/cazy_readmap/SAMPLE_hits.tsv [more_hits.tsv ...]

suppressPackageStartupMessages(library(data.table))
options(width = 150)

PIDENT_MIN <- 50
EVALUE_MAX <- 1e-9
LABEL_PAT  <- "^(GH|GT|PL|CE|AA|CBM)[0-9]+(_[0-9]+)?$"
BIN_PAT    <- "^(GH|GT|PL|CE|AA|CBM)0$"
EC_PAT     <- "^[0-9]+\\.[0-9-]+"

base_family <- function(x) sub("_[0-9]+$", "", x)

check_hits <- function(hits_file) {
  cat("============================================================\n")
  cat("File:", hits_file, "\n")
  cat("============================================================\n")

  dt <- fread(hits_file, header = FALSE, sep = "\t", quote = "", fill = Inf,
              select = c(1L, 2L, 3L, 11L, 12L),
              col.names = c("qseqid", "sseqid", "pident", "evalue", "bitscore"))
  dt <- dt[pident >= PIDENT_MIN & evalue <= EVALUE_MAX]
  cat("Filtered hits:", nrow(dt), "\n")
  if (nrow(dt) == 0) return(invisible(NULL))

  # Parse each distinct header once
  u      <- unique(dt$sseqid)
  fields <- strsplit(u, "|", fixed = TRUE)
  labs   <- lapply(fields, function(p) p[grepl(LABEL_PAT, p)])
  lu <- data.table(
    sseqid   = u,
    field1   = sapply(fields, `[`, 1L),
    n_bases  = sapply(labs, function(x) length(unique(base_family(x)))),
    first    = sapply(labs, function(x) if (length(x)) x[1] else NA_character_),
    last     = sapply(labs, function(x) if (length(x)) x[length(x)] else NA_character_),
    has_cat  = sapply(labs, function(x) any(!grepl("^CBM", x))),
    has_ec   = sapply(fields, function(p) any(grepl(EC_PAT, p)))
  )
  lu[, genbank := grepl(".", field1, fixed = TRUE)]
  # The label 08_integrate.R keeps: last for GenBank headers, first otherwise
  lu[, used := fifelse(genbank, last, first)]
  dt <- lu[dt, on = "sseqid"]

  # ── 1. Multi-family headers ────────────────────────────────────────────────
  cat("\n--- 1. Hits by number of distinct families in the header ---\n")
  print(dt[, .N, by = n_bases][order(n_bases)][, pct := round(100 * N / sum(N), 2)][])
  cat("\nMulti-family hits by header format:\n")
  print(dt[n_bases > 1, .N, by = .(format = fifelse(genbank, "GenBank", "JGI/other"))])
  n_cbm <- dt[grepl("^CBM", used) & has_cat == TRUE, .N]
  cat(sprintf("\nHits kept as CBM although the header has a catalytic family: %d (%.2f%%)\n",
              n_cbm, 100 * n_cbm / nrow(dt)))
  cat("Example multi-family headers:\n")
  print(head(lu[n_bases > 1, sseqid], 10))

  # ── 2. Unassigned bins ─────────────────────────────────────────────────────
  cat("\n--- 2. Unassigned family 0 bins vs. classified families ---\n")
  dt[, group := fifelse(grepl(BIN_PAT, used), used, "classified")]
  print(dt[!is.na(used), .(
    hits           = .N,
    pct_hits       = round(100 * .N / dt[!is.na(used), .N], 2),
    med_pident     = round(median(pident), 1),
    med_bitscore   = round(median(bitscore), 1),
    pct_with_ec    = round(100 * mean(has_ec), 1),
    pct_other_fam  = round(100 * mean(n_bases > 1), 1),
    distinct_refs  = uniqueN(sseqid)
  ), by = group][order(-hits)])
  cat("\nMost-hit references in each bin:\n")
  print(dt[group != "classified", .N, by = .(group, sseqid)][
    order(group, -N)][, head(.SD, 3), by = group])

  # ── 3. Unrecognized header formats ─────────────────────────────────────────
  # No version dot sends these down the JGI branch, where field 1 is taken as the
  # family, so an unversioned accession becomes a false family name.
  cat("\n--- 3. Headers with an unrecognized first field ---\n")
  odd <- dt[genbank == FALSE & !grepl(LABEL_PAT, field1), .N, by = sseqid][order(-N)]
  cat("Hits:", sum(odd$N), " distinct headers:", nrow(odd), "\n")
  print(head(odd, 20))

  # ── 4. Paired-read double counts ───────────────────────────────────────────
  cat("\n--- 4. Read-pair counts ---\n")
  counted <- dt[!is.na(used)]
  counted[, pair_id := sub("[/. ][12]$", "", qseqid)]
  raw    <- uniqueN(counted[, .(pair_id, used)])
  merged <- uniqueN(counted[, .(pair_id, base_family(used))])
  pairs  <- uniqueN(counted$pair_id)
  cat("Dedup on labels as written      :", raw, "\n")
  cat("Merge to family, then dedup     :", merged, "\n")
  cat(sprintf("Double counts removed by merging: %d (%.2f%%)\n",
              raw - merged, 100 * (raw - merged) / merged))
  cat(sprintf("Distinct read pairs             : %d\n", pairs))
  cat(sprintf("Pairs still in >1 family        : %d (%.2f%% of pairs)\n\n",
              merged - pairs, 100 * (merged - pairs) / pairs))
  invisible(NULL)
}

files <- commandArgs(trailingOnly = TRUE)
if (length(files) == 0) stop("Usage: Rscript check_hits.R SAMPLE_hits.tsv [more_hits.tsv ...]", call. = FALSE)
for (f in files) check_hits(f)
