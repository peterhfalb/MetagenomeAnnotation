#!/usr/bin/env Rscript
# check_family_refs.R
# Read-only check of whether fungal:bacterial ratios for a CAZy family can be
# trusted, using the database and one sample's hits:
#   1. How many reference sequences each family has, by kingdom. A family with no
#      bacterial references can't produce bacterial counts, whatever is in the soil.
#   2. For selected families, which fungal references the reads hit, and at what
#      identity. A few heavily hit references can point to contamination in a
#      genome or to reads from another organism.
#
# Needs check_cazy_db.sh output (headers.txt, taxid_lineages.tsv).
# Writes family_reference_counts.tsv to the check directory.
#
# Usage:
#   Rscript check_family_refs.R SAMPLE_hits.tsv [check_dir] [taxonmap]

suppressPackageStartupMessages(library(data.table))
options(width = 150, datatable.print.nrows = 200)

args      <- commandArgs(trailingOnly = TRUE)
hits_file <- args[1]
check_dir <- if (length(args) >= 2) args[2] else file.path(Sys.getenv("HOME"), "cazy_db_check")
taxonmap  <- if (length(args) >= 3) args[3] else
  "/projects/standard/kennedyp/shared/databases/metaG_annotation/cazy_readmap/cazy_taxonmap.tsv"
if (is.na(hits_file)) stop("Usage: Rscript check_family_refs.R SAMPLE_hits.tsv [check_dir] [taxonmap]", call. = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
jgi_file   <- file.path(dirname(sub("^--file=", "", script_arg[1])), "jgi_fungal_genomes.tsv")

LABEL_PAT <- "^(GH|GT|PL|CE|AA|CBM)[0-9]+(_[0-9]+)?$"
FOCUS     <- c("AA1", "AA2", "AA3", "AA4", "AA5", "AA6", "AA9", "GH23", "GH24", "GH25", "CE5", "GH93")
base_family <- function(x) sub("_[0-9]+$", "", x)

tm  <- fread(taxonmap, header = TRUE, col.names = c("accession", "taxid"),
             colClasses = c("character", "integer"), key = "accession")
lin <- fread(file.path(check_dir, "taxid_lineages.tsv"), sep = "\t", quote = "",
             select = c("taxid", "group"), colClasses = list(integer = "taxid", character = "group"))
jgi <- fread(jgi_file, sep = "\t", quote = "", select = "source_file")$source_file

# Kingdom group of each reference, by the same rules as 08_integrate.R
ref_group <- function(sseqid) {
  field1 <- sub("\\|.*$", "", sseqid)
  is_jgi <- grepl(LABEL_PAT, field1)
  out <- rep(NA_character_, length(sseqid))
  out[is_jgi] <- fifelse(sub("^.*\\|", "", sseqid[is_jgi]) %in% jgi, "JGI_Fungi", "JGI_unverified")
  taxid <- tm$taxid[match(field1[!is_jgi], tm$accession)]
  g <- lin$group[match(taxid, lin$taxid)]
  g[is.na(taxid)] <- "no_taxid"
  g[is.na(g)] <- "unresolved"
  out[!is_jgi] <- paste0("NCBI_", g)
  out
}

# ── 1. Reference sequences per family, by kingdom ────────────────────────────
headers <- readLines(file.path(check_dir, "headers.txt"), warn = FALSE)
cat("Database headers:", length(headers), "\n")
groups <- ref_group(headers)
fields <- strsplit(headers, "|", fixed = TRUE)
long <- data.table(seq = rep(seq_along(fields), lengths(fields)), label = unlist(fields, use.names = FALSE))
long <- unique(long[grepl(LABEL_PAT, label), .(seq, family = base_family(label))])
long[, group := groups[seq]]
refs <- dcast(long, family ~ group, value.var = "seq", fun.aggregate = length)
fwrite(refs, file.path(check_dir, "family_reference_counts.tsv"), sep = "\t")

cols <- setdiff(names(refs), "family")
fungal   <- intersect(cols, c("JGI_Fungi", "NCBI_Fungi"))
refs[, fungal_refs := rowSums(.SD), .SDcols = fungal]
refs[, bacterial_refs := if ("NCBI_Bacteria" %in% cols) NCBI_Bacteria else 0L]
refs[, other_refs := rowSums(.SD), .SDcols = setdiff(cols, c(fungal, "NCBI_Bacteria"))]

cat("\n--- 1. Reference sequences for selected families ---\n")
print(refs[family %in% FOCUS, .(family, fungal_refs, bacterial_refs, other_refs)][order(family)])
cat("\nFamilies with fungal references and no bacterial references:",
    refs[fungal_refs > 0 & bacterial_refs == 0, .N], "of", nrow(refs), "\n")
print(refs[fungal_refs > 0 & bacterial_refs == 0, .(family, fungal_refs, other_refs)][order(-fungal_refs)][1:min(.N, 40)])
cat("\nFamilies with no fungal references:", refs[fungal_refs == 0, .N], "\n")

# ── 2. What reads hit in one sample ──────────────────────────────────────────
dt <- fread(hits_file, header = FALSE, sep = "\t", quote = "", fill = Inf,
            select = c(2L, 3L, 11L), col.names = c("sseqid", "pident", "evalue"))
dt <- dt[pident >= 50 & evalue <= 1e-9]
dt[, pident := as.numeric(pident)]
u  <- unique(dt$sseqid)
pick <- vapply(strsplit(u, "|", fixed = TRUE), function(p) {
  labs <- p[grepl(LABEL_PAT, p)]
  if (length(labs) == 0) return(NA_character_)
  is_bin <- grepl("^[A-Z]+0$", labs)
  is_cbm <- grepl("^CBM", labs)
  c(labs[!is_bin & !is_cbm], labs[!is_bin & is_cbm], labs[is_bin])[1]
}, character(1))
ug <- ref_group(u)
dt[, `:=`(family = base_family(pick[match(sseqid, u)]), group = ug[match(sseqid, u)])]
dt <- dt[!is.na(family)]

cat("\n--- 2. Identity of hits by reference group (", basename(hits_file), ") ---\n", sep = "")
print(dt[, .(hits = .N, q25 = quantile(pident, 0.25), median = median(pident), q75 = quantile(pident, 0.75),
             pct_ge70 = round(100 * mean(pident >= 70), 1), pct_ge90 = round(100 * mean(pident >= 90), 1)),
         by = group][order(-hits)])

cat("\nIdentity of fungal hits in selected families:\n")
print(dt[family %in% FOCUS & group %in% c("JGI_Fungi", "NCBI_Fungi"),
         .(hits = .N, median = median(pident), pct_ge70 = round(100 * mean(pident >= 70), 1),
           pct_ge90 = round(100 * mean(pident >= 90), 1)), by = family][order(family)])

cat("\nMost-hit fungal references in selected families (genome file or accession):\n")
fungal_hits <- dt[family %in% FOCUS & group %in% c("JGI_Fungi", "NCBI_Fungi")]
fungal_hits[, ref := fifelse(group == "JGI_Fungi", sub("_(GeneCatalog|FilteredModels|filtered).*$", "", sub("^.*\\|", "", sseqid)), sseqid)]
top <- fungal_hits[, .(hits = .N), by = .(family, ref)][order(family, -hits)]
top[, pct_of_family := round(100 * hits / sum(hits), 1), by = family]
print(top[, head(.SD, 4), by = family])
