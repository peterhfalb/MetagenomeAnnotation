#!/usr/bin/env Rscript
# update_summary_stats.R
# Add or refresh stand_type in an existing summary_stats.tsv without re-running
# the full 08_integrate.R pipeline. Reads the TSV, joins stakes.csv on the
# extracted stake ID (GEO-{row}-{stake}_... -> row-stake), and overwrites the file.
#
# Usage:
#   Rscript update_summary_stats.R \
#     --summary   /path/to/integrated/summary_stats.tsv \
#     --metadata  /path/to/stakes.csv

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--summary",  type = "character", default = NULL,
              help = "Path to summary_stats.tsv (required)"),
  make_option("--metadata", type = "character", default = NULL,
              help = "Path to stakes.csv with columns: stake, stand_type (required)")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$summary) || is.null(opt$metadata))
  stop("--summary and --metadata are required.", call. = FALSE)
if (!file.exists(opt$summary))  stop("summary_stats.tsv not found: ",  opt$summary,  call. = FALSE)
if (!file.exists(opt$metadata)) stop("stakes.csv not found: ",          opt$metadata, call. = FALSE)

dt   <- fread(opt$summary,  colClasses = "character")
meta <- fread(opt$metadata, colClasses = "character")

if (!all(c("stake", "stand_type") %in% names(meta)))
  stop("stakes.csv must have columns 'stake' and 'stand_type'.", call. = FALSE)

# Drop any existing stake_id / stand_type columns so we can re-join cleanly
dt[, c("stake_id", "stand_type") := NULL]

# Extract stake (e.g. GEO-00-06_S1 -> 00-06, GEO-03-13_Cleaned_S8 -> 03-13)
dt[, stake_id := sub("^GEO-([0-9]+-[0-9]+).*$", "\\1", sample)]

dt <- merge(dt, meta[, .(stake_id = stake, stand_type)],
            by = "stake_id", all.x = TRUE, sort = FALSE)

setcolorder(dt, c("sample", "stake_id", "stand_type"))

n_matched <- sum(!is.na(dt$stand_type))
cat("stand_type joined for", n_matched, "of", nrow(dt), "samples\n")
if (n_matched < nrow(dt))
  cat("WARNING:", nrow(dt) - n_matched, "samples had no match in stakes.csv\n")

fwrite(dt, opt$summary, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote:", opt$summary, "\n")
