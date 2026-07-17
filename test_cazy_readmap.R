#!/usr/bin/env Rscript
# test_cazy_readmap.R
# Quick end-to-end test of the CAZy readmap taxonomy pipeline on MSI.
# Run after the test job (submit_cazy_readmap.sh --test 1) completes.
#
# Usage:
#   module load R
#   Rscript test_cazy_readmap.R

suppressPackageStartupMessages({
  library(data.table)
  library(taxonomizr)
})

DB_DIR  <- "/projects/standard/kennedyp/shared/databases/metaG_annotation"
ANN_DIR <- "/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation"

cat("============================================================\n")
cat("CAZy readmap taxonomy pipeline — test run\n")
cat("============================================================\n\n")

# --- 1. Locate hits file ---
sample_list <- file.path(ANN_DIR, "sample_list.txt")
sample <- readLines(sample_list, n = 1)
hits_file <- file.path(ANN_DIR, "cazy_readmap", paste0(sample, "_hits.tsv"))

cat("Sample :", sample, "\n")

# Accept partial R1 file if the merged file isn't ready yet
r1_file <- file.path(ANN_DIR, "cazy_readmap", paste0(sample, "_R1_hits.tsv"))
if (!file.exists(hits_file) && file.exists(r1_file)) {
  cat("Note   : merged hits not ready — using partial R1 file for testing\n")
  hits_file <- r1_file
}
cat("Hits   :", hits_file, "\n")
if (!file.exists(hits_file)) stop("Hits file not found — did the SLURM job finish?")
cat("Size   :", format(file.info(hits_file)$size, big.mark = ","), "bytes\n\n")

# --- 2. Check supporting files ---
taxonmap_file <- file.path(DB_DIR, "cazy_readmap", "cazy_taxonmap.tsv")
taxa_sql      <- file.path(DB_DIR, "taxonomy", "taxonomy.sql")

cat("Taxonmap  :", ifelse(file.exists(taxonmap_file), "[OK]", "[MISSING]"), taxonmap_file, "\n")
cat("Taxonomy  :", ifelse(file.exists(taxa_sql),      "[OK]", "[MISSING — will be built on first full 08_integrate.R run]"), "\n\n")

# --- 3. Parse hits ---
cat("Parsing hits...\n")
dt <- fread(hits_file, header = FALSE, sep = "\t", quote = "", fill = Inf,
            select = c(1L, 2L, 3L, 11L, 12L),
            col.names = c("qseqid", "sseqid", "pident", "evalue", "bitscore"))

cat("  Total hits (raw)          :", formatC(nrow(dt), format = "d", big.mark = ","), "\n")
dt <- dt[pident >= 50 & evalue <= 1e-9]
cat("  After pident>=50 e<=1e-9  :", formatC(nrow(dt), format = "d", big.mark = ","), "\n")

FAMILY_PAT <- "([A-Za-z]+[0-9]+(?:_[0-9]+)?)"
dt[, field1     := sub("\\|.*$", "", sseqid)]
dt[, is_genbank := grepl("\\.", field1)]

dt[, family := {
  f <- rep(NA_character_, .N)
  gb  <- is_genbank
  f[gb]  <- sub(paste0(".*\\|", FAMILY_PAT, "(?:\\|[^A-Za-z].*)?$"), "\\1", sseqid[gb])
  f[gb & f == sseqid] <- NA_character_
  jgi <- !is_genbank
  f[jgi] <- sub(paste0("^", FAMILY_PAT, ".*$"), "\\1", sseqid[jgi])
  f[jgi & !grepl("^[A-Za-z]+[0-9]", f)] <- NA_character_
  f
}]
dt <- dt[!is.na(family)]
cat("  After family extraction    :", formatC(nrow(dt), format = "d", big.mark = ","), "\n")

dt[is_genbank == TRUE,  accession := field1]
dt[is_genbank == FALSE, taxid := {
  parts <- strsplit(sseqid[is_genbank == FALSE], "\\|")
  as.integer(sapply(parts, function(p) {
    num <- p[grepl("^[0-9]+$", p)]
    if (length(num) > 0) tail(num, 1L) else NA_character_
  }))
}]
dt[, c("field1", "is_genbank") := NULL]

dt[, pair_id := sub("[/. ][12]$", "", qseqid)]
dt <- dt[dt[, .I[which.max(bitscore)], by = .(pair_id, family)]$V1]
cat("  After paired-read dedup   :", formatC(nrow(dt), format = "d", big.mark = ","), "\n\n")

# --- 4. Taxonmap join ---
cat("Joining taxonmap...\n")
tm <- fread(taxonmap_file, header = TRUE,
            col.names = c("accession", "taxid_map"),
            colClasses = c("character", "integer"), key = "accession")
needs <- !is.na(dt$accession)
joined <- tm[dt[needs, .(accession)], on = "accession"]
dt[needs, taxid := joined$taxid_map]
dt[, accession := NULL]

cat("  Taxids assigned            :", formatC(sum(!is.na(dt$taxid)), format = "d", big.mark = ","),
    sprintf("/ %s (%.1f%%)\n",
            formatC(nrow(dt), format = "d", big.mark = ","),
            100 * sum(!is.na(dt$taxid)) / nrow(dt)))
cat("  No taxid (unclassified)   :", formatC(sum(is.na(dt$taxid)), format = "d", big.mark = ","), "\n\n")

# --- 5. Kingdom expansion ---
if (!file.exists(taxa_sql)) {
  cat("taxonomy.sql not yet built — skipping kingdom expansion.\n")
  cat("It will be built automatically on the first full run of 08_integrate.R.\n\n")
} else {
  cat("Expanding taxids to kingdom via taxonomizr...\n")
  all_taxids <- unique(dt$taxid[!is.na(dt$taxid) & dt$taxid > 0])
  cat("  Unique taxids:", length(all_taxids), "\n")

  taxon_mat <- getTaxonomy(all_taxids, taxa_sql)

  cat("  Columns returned by getTaxonomy():", paste(colnames(taxon_mat), collapse = ", "), "\n")

  # Pick the superkingdom column — name varies across taxonomizr versions
  sk_col <- intersect(c("superkingdom", "domain", "Superkingdom", "super kingdom"), colnames(taxon_mat))[1]
  kk_col <- intersect(c("kingdom",      "Kingdom"),                                colnames(taxon_mat))[1]
  if (is.na(sk_col)) stop("Cannot find superkingdom column in getTaxonomy() output — columns: ",
                          paste(colnames(taxon_mat), collapse = ", "))

  FUNGAL_PHYLA <- c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mucoromycota",
                    "Glomeromycota", "Blastocladiomycota", "Neocallimastigomycota",
                    "Zoopagomycota", "Mortierellomycota", "Microsporidia")

  kingdom_lut <- data.table(
    taxid  = as.integer(rownames(taxon_mat)),
    domain = taxon_mat[, sk_col],
    phylum = if ("phylum" %in% colnames(taxon_mat)) taxon_mat[, "phylum"] else NA_character_
  )
  # Assign kingdom: use phylum to distinguish Fungi from other Eukaryota
  kingdom_lut[, kingdom := domain]
  if (!is.na(kk_col)) {
    kingdom_lut[!is.na(taxon_mat[, kk_col]) & taxon_mat[, kk_col] == "Fungi",
                kingdom := "Fungi"]
  } else {
    kingdom_lut[phylum %in% FUNGAL_PHYLA, kingdom := "Fungi"]
  }
  kingdom_lut[is.na(kingdom), kingdom := "unclassified"]

  dt[kingdom_lut, kingdom := i.kingdom, on = "taxid"]
  dt[is.na(kingdom), kingdom := "unclassified"]

  cat("\nKingdom breakdown:\n")
  ks <- dt[, .N, by = kingdom][order(-N)]
  ks[, pct := round(100 * N / sum(N), 1)]
  print(ks)

  n_f   <- ks[kingdom == "Fungi", N]
  pct_f <- ks[kingdom == "Fungi", pct]
  cat(sprintf("\n>>> Fungi: %s reads (%.1f%% of %s filtered, deduplicated hits)\n\n",
              formatC(n_f, format = "d", big.mark = ","),
              pct_f,
              formatC(nrow(dt), format = "d", big.mark = ",")))

  cat("Top 10 Fungi families:\n")
  fungi_fam <- dt[kingdom == "Fungi", .N, by = family][order(-N)]
  if (nrow(fungi_fam) == 0) cat("  (no Fungi hits — check kingdom assignment)\n")
  else print(fungi_fam[1:min(10, .N)])
}

cat("\n============================================================\n")
cat("Test complete.\n")
cat("============================================================\n")
