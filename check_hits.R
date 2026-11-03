#!/usr/bin/env Rscript
# check_hits.R
# Diagnostics on raw CAZy readmap hits (01c_cazy_readmap.sh output). Read-only;
# applies the same pident/e-value filters as 08_integrate.R, then reports:
#   1. Multi-family headers, where the label kept depends on label order
#   2. Unassigned "family 0" bins (GH0, AA0, ...)
#   3. Headers whose first field is neither a CAZy label nor a versioned accession
#   4. Paired-read double counts removed by merging to family before dedup
#   5. Kingdoms assigned to JGI-format references, whose numeric header field
#      08_integrate.R treats as an NCBI taxid. Writes {SAMPLE}_jgi_codes.tsv to
#      the working directory; everything else only prints.
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

DB_DIR        <- "/projects/standard/kennedyp/shared/databases/metaG_annotation"
TAXONMAP_FILE <- file.path(DB_DIR, "cazy_readmap", "cazy_taxonmap.tsv")
TAXA_SQL      <- file.path(DB_DIR, "taxonomy", "taxonomy.sql")
FUNGAL_PHYLA  <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mucoromycota",
                   "Glomeromycota", "Blastocladiomycota", "Neocallimastigomycota",
                   "Zoopagomycota", "Mortierellomycota", "Microsporidia")

base_family <- function(x) sub("_[0-9]+$", "", x)

# Same taxid -> kingdom logic as 08_integrate.R
kingdom_lookup <- function(taxids) {
  if (!requireNamespace("taxonomizr", quietly = TRUE)) return(NULL)
  mat <- taxonomizr::getTaxonomy(taxids, TAXA_SQL)
  sk  <- intersect(c("superkingdom", "domain", "Superkingdom", "super kingdom"), colnames(mat))[1]
  kk  <- intersect(c("kingdom", "Kingdom"), colnames(mat))[1]
  if (is.na(sk)) return(NULL)
  out <- data.table(taxid = as.integer(rownames(mat)), kingdom = mat[, sk])
  if (!is.na(kk)) {
    out[!is.na(mat[, kk]) & mat[, kk] == "Fungi", kingdom := "Fungi"]
  } else if ("phylum" %in% colnames(mat)) {
    out[mat[, "phylum"] %in% FUNGAL_PHYLA, kingdom := "Fungi"]
  }
  out[is.na(kingdom), kingdom := "unclassified"]
  out
}

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
    has_ec   = sapply(fields, function(p) any(grepl(EC_PAT, p))),
    # Last all-numeric field, as 08_integrate.R reads it for JGI headers
    num_id   = sapply(fields, function(p) {
      num <- p[grepl("^[0-9]+$", p)]
      if (length(num)) num[length(num)] else NA_character_
    }),
    last_fld = sapply(fields, function(p) p[length(p)])
  )
  lu[, genbank := grepl(".", field1, fixed = TRUE)]
  lu[, jgi := !genbank & grepl(LABEL_PAT, field1)]
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

  # ── 5. Kingdoms assigned to JGI references ─────────────────────────────────
  cat("--- 5. Kingdom assignment by reference source ---\n")
  counted[, family := base_family(used)]
  k <- counted[counted[, .I[which.max(bitscore)], by = .(pair_id, family)]$V1]
  k[, source := fifelse(jgi, "JGI", fifelse(genbank, "GenBank", "other"))]

  codes <- k[jgi == TRUE, .(hits = .N, file_field = last_fld[1]),
             by = .(code = sub("_.*$", "", last_fld))][order(-hits)]
  codes_file <- paste0(sub("_hits\\.tsv$", "", basename(hits_file)), "_jgi_codes.tsv")
  fwrite(codes, codes_file, sep = "\t")
  cat("JGI counts:", sum(codes$hits), " distinct genome codes:", nrow(codes),
      " written to:", codes_file, "\n")
  cat("Most-hit JGI genomes:\n")
  print(head(codes[, .(code, hits)], 15))
  cat("Share of JGI counts from GeneCatalog files:",
      round(100 * k[jgi == TRUE, mean(grepl("GeneCatalog", last_fld))], 1), "%\n")

  if (!file.exists(TAXONMAP_FILE) || !file.exists(TAXA_SQL)) {
    cat("\n[SKIP] taxonmap or taxonomy.sql not found; kingdom comparison not run.\n\n")
    return(invisible(NULL))
  }
  tm <- fread(TAXONMAP_FILE, header = TRUE, col.names = c("accession", "taxid"),
              colClasses = c("character", "integer"), key = "accession")
  k[, taxid := NA_integer_]
  k[jgi == TRUE, taxid := as.integer(num_id)]
  k[genbank == TRUE, taxid := tm[.SD, on = c(accession = "field1"), x.taxid]]

  ids <- unique(k$taxid[!is.na(k$taxid) & k$taxid > 0])
  kl  <- kingdom_lookup(ids)
  if (is.null(kl)) {
    cat("\n[SKIP] Could not resolve kingdoms; is taxonomizr installed?\n\n")
    return(invisible(NULL))
  }
  k[kl, kingdom := i.kingdom, on = "taxid"]
  k[is.na(kingdom), kingdom := "unclassified"]

  cat("\nCounts by reference source and assigned kingdom:\n")
  print(dcast(k, source ~ kingdom, value.var = "pair_id", fun.aggregate = length))

  # Upper bound: treats every JGI reference as fungal, which needs confirming
  # against the genome codes above before it is used for anything.
  fb <- function(d) {
    f <- d[kingdom == "Fungi", .N]; b <- d[kingdom == "Bacteria", .N]
    data.table(fungi = f, bacteria = b, log_ratio = round(log((f + 1) / (b + 1)), 3))
  }
  all_jgi_fungal <- copy(k)[jgi == TRUE, kingdom := "Fungi"]
  cat("\nF:B as parsed now vs. if every JGI reference were fungal:\n")
  print(rbind(cbind(scenario = "as parsed", fb(k)),
              cbind(scenario = "JGI = Fungi", fb(all_jgi_fungal))))
  focus <- c("AA0", "AA1", "AA3", "AA5", "AA9", "GH5", "GH43", "GH16")
  cat("\nSame comparison for selected families:\n")
  print(rbindlist(lapply(focus, function(fam) {
    cbind(family = fam,
          as_parsed = fb(k[family == fam]),
          jgi_fungi = fb(all_jgi_fungal[family == fam]))
  })))
  cat("\nFamilies with the most JGI counts assigned to Bacteria:\n")
  print(head(k[jgi == TRUE & kingdom == "Bacteria", .N, by = family][order(-N)], 15))
  cat("\n")
  invisible(NULL)
}

files <- commandArgs(trailingOnly = TRUE)
if (length(files) == 0) stop("Usage: Rscript check_hits.R SAMPLE_hits.tsv [more_hits.tsv ...]", call. = FALSE)
for (f in files) {
  if (!file.exists(f)) {
    cat("[SKIP] File not found (check the full path):", f, "\n\n")
    next
  }
  check_hits(f)
}
