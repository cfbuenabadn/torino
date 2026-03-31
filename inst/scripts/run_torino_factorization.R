#!/usr/bin/env Rscript
# Run Torino factorization on an input count matrix.
#
# Usage:
#   Rscript run_torino_factorization.R <input_csv> <min_counts> <K> <out_rds> [merge_n]
#
# Positional arguments:
#   input_csv   CSV file with a Sample_ID column (rows = samples, remaining
#               columns = genomic positions).
#   min_counts  Minimum total read count per sample. Samples below this are
#               removed before factorization.
#   K           Number of factors.
#   out_rds     Output RDS file path.
#   merge_n     (optional) Number of consecutive position columns to merge to
#               reduce RAM. Use 'auto' (default) for automatic size-based
#               selection, or an integer such as 1 (no merging), 3, 9, or 30.
#
# Optional sample filtering (environment variable):
#   TORINO_FILTER_SAMPLES   Path to a plain-text file with one sample ID per
#                           line. Only samples present in this list are kept.
#
# Output RDS structure:
#   $gene        Gene name derived from the input filename.
#   $torino      List returned by torino_factorization() (EF, EF_smooth, EL,
#                elbo, coords, samples), or NULL if factorization failed.
#   $K           Number of factors requested.
#   $coords_all  All position column names from the original (unmerged) matrix.

suppressPackageStartupMessages({
    library(torino)
    library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    cat(
        "Usage: Rscript run_torino_factorization.R",
        "<input_csv> <min_counts> <K> <out_rds> [merge_n]\n"
    )
    quit(status = 1)
}

input_matrix <- args[1]
min_counts   <- as.integer(args[2])
K            <- as.integer(args[3])
out_file     <- args[4]
merge_n_arg  <- if (length(args) >= 5) args[5] else "auto"

gene_name <- basename(input_matrix) %>% tools::file_path_sans_ext()

message(sprintf("Loading %s ...", input_matrix))
counts <- read_csv(input_matrix, show_col_types = FALSE) %>%
    column_to_rownames(var = "Sample_ID")
# Strip version suffixes from sample IDs (e.g. SAMPLE.1 -> SAMPLE)
rownames(counts) <- sub("\\.[^.]*$", "", rownames(counts))

coords_all   <- colnames(counts)
nbases_total <- ncol(counts)

# ── Column merging ────────────────────────────────────────────────────────────
if (merge_n_arg == "auto") {
    if (nbases_total > 1000000) {
        n <- 30
        message(sprintf("Matrix has %d bases (>1M); merging groups of 30.", nbases_total))
    } else if (nbases_total > 60000) {
        n <- 9
        message(sprintf("Matrix has %d bases (>60k); merging groups of 9.", nbases_total))
    } else if (nbases_total > 15000) {
        n <- 3
        message(sprintf("Matrix has %d bases (>15k); merging groups of 3.", nbases_total))
    } else {
        n <- 1
        message(sprintf("Matrix has %d bases; no column merging needed.", nbases_total))
    }
} else {
    n <- as.integer(merge_n_arg)
    if (n > 1) message(sprintf("Merging groups of %d columns (user-specified).", n))
}

# ── Optional sample filtering ─────────────────────────────────────────────────
filter_file <- Sys.getenv("TORINO_FILTER_SAMPLES", unset = "")
if (nzchar(filter_file)) {
    select_samples <- readLines(filter_file)
    keep <- intersect(select_samples, rownames(counts))
    if (length(keep) == 0) stop("No samples from TORINO_FILTER_SAMPLES found in the count matrix.")
    counts <- counts[keep, , drop = FALSE]
    message(sprintf("Filtered to %d / %d samples using %s",
                    length(keep), length(select_samples), filter_file))
}

# ── Preprocessing ─────────────────────────────────────────────────────────────
message(sprintf("Preparing counts (min_counts=%d, n=%d) ...", min_counts, n))
counts <- prepare_counts(counts, min_counts, n = n)
message(sprintf("Matrix after preparation: %d samples x %d positions", nrow(counts), ncol(counts)))

# ── Factorization ─────────────────────────────────────────────────────────────
message(sprintf("Running torino_factorization with K=%d ...", K))
fit <- torino_factorization(counts, K = K, seed = 1)

if (is.null(fit)) {
    warning("torino_factorization returned NULL (factorization failed).")
} else {
    message("Factorization complete.")
}

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(
    list(
        gene       = gene_name,
        torino     = fit,
        K          = K,
        coords_all = coords_all
    ),
    file = out_file
)
message(sprintf("Results saved to %s", out_file))
