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
nodes_file    <- file.path(DB_DIR, "taxonomy", "nodes.dmp")

cat("Taxonmap  :", ifelse(file.exists(taxonmap_file), "[OK]", "[MISSING]"), taxonmap_file, "\n")
cat("Taxonomy  :", ifelse(file.exists(nodes_file),    "[OK]", "[MISSING]"), nodes_file, "\n\n")

# --- 3. Parse hits ---
cat("Parsing hits...\n")
dt <- fread(hits_file, header = FALSE, sep = "\t", quote = "", fill = Inf,
            select = c(1L, 2L, 3L, 11L, 12L),
            col.names = c("qseqid", "sseqid", "pident", "evalue", "bitscore"))

cat("  Total hits (raw)          :", formatC(nrow(dt), format = "d", big.mark = ","), "\n")
dt <- dt[pident >= 50 & evalue <= 1e-9]
cat("  After pident>=50 e<=1e-9  :", formatC(nrow(dt), format = "d", big.mark = ","), "\n")

LABEL_PAT <- "^(GH|GT|PL|CE|AA|CBM)[0-9]+(_[0-9]+)?$"
dt[, field1     := sub("\\|.*$", "", sseqid)]
dt[, is_genbank := !grepl(LABEL_PAT, field1)]

sseqids <- unique(dt$sseqid)
labels  <- vapply(strsplit(sseqids, "|", fixed = TRUE), function(p) {
  labs <- p[grepl(LABEL_PAT, p)]
  if (length(labs) == 0) return(NA_character_)
  is_bin <- grepl("^[A-Z]+0$", labs)
  is_cbm <- grepl("^CBM", labs)
  c(labs[!is_bin & !is_cbm], labs[!is_bin & is_cbm], labs[is_bin])[1]
}, character(1))
dt[, family := labels[match(sseqid, sseqids)]]
dt <- dt[!is.na(family)]
cat("  After family extraction    :", formatC(nrow(dt), format = "d", big.mark = ","), "\n")

# JGI numbers are protein IDs, not taxids: the genome file identifies the organism
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
jgi_genomes <- fread(file.path(dirname(sub("^--file=", "", script_arg[1])), "jgi_fungal_genomes.tsv"),
                     sep = "\t", quote = "", select = c("source_file", "ncbi_taxid"),
                     colClasses = list(character = "source_file", integer = "ncbi_taxid"))
dt[, taxid := NA_integer_]
dt[is_genbank == TRUE,  accession := field1]
dt[is_genbank == FALSE, taxid := jgi_genomes$ncbi_taxid[
  match(sub("^.*\\|", "", sseqid), jgi_genomes$source_file)]]
dt[, source := fifelse(is_genbank, "NCBI", "JGI")]
dt[, c("field1", "is_genbank") := NULL]

dt[, family := sub("_[0-9]+$", "", family)]
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
if (!file.exists(nodes_file)) {
  cat("nodes.dmp not found — skipping kingdom assignment.\n\n")
} else {
  cat("Assigning kingdoms from NCBI lineages...\n")
  all_taxids <- unique(dt$taxid[!is.na(dt$taxid) & dt$taxid > 0])
  cat("  Unique taxids:", length(all_taxids), "\n")

  nodes <- fread(nodes_file, sep = "\t", header = FALSE, quote = "",
                 select = c(1L, 3L), col.names = c("taxid", "parent"))
  tops <- c("4751" = "Fungi", "2" = "Bacteria", "2157" = "Archaea",
            "10239" = "Viruses", "2759" = "Eukaryota")
  kingdoms <- rep("unclassified", length(all_taxids))
  current  <- all_taxids
  for (step in seq_len(100L)) {
    active <- which(!is.na(current))
    if (length(active) == 0) break
    hit   <- unname(tops[as.character(current[active])])
    found <- !is.na(hit)
    kingdoms[active[found]] <- hit[found]
    current[active[found]]  <- NA_integer_
    rest   <- active[!found]
    parent <- nodes$parent[match(current[rest], nodes$taxid)]
    parent[parent == current[rest]] <- NA_integer_
    current[rest] <- parent
  }
  kingdom_lut <- data.table(taxid = all_taxids, kingdom = kingdoms)

  dt[kingdom_lut, kingdom := i.kingdom, on = "taxid"]
  dt[is.na(kingdom), kingdom := "unclassified"]

  cat("\nCounts by reference source and kingdom:\n")
  print(dt[, .N, by = .(source, kingdom)][order(source, -N)])

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
