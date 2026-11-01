#' Pick a column-merging factor from the matrix width
#'
#' Implements the size-based rule used by the Torino pipeline: wide loci are
#' binned more aggressively so that the factorization fits in RAM.
#'
#' @param nbases Integer. Number of position columns in the count matrix.
#' @return Integer merging factor: 30 for >1M bases, 9 for >60k, 3 for >15k,
#'   otherwise 1 (no merging).
#' @export
choose_merge_n <- function(nbases) {
    if (nbases > 1000000) {
        30
    } else if (nbases > 60000) {
        9
    } else if (nbases > 15000) {
        3
    } else {
        1
    }
}


#' Read a Torino count CSV
#'
#' Reads a (possibly gzipped) CSV whose first column is \code{Sample_ID} and
#' whose remaining columns are genomic positions, and moves \code{Sample_ID}
#' into the row names. Raises \code{VROOM_CONNECTION_SIZE} for the duration of
#' the read, since these headers are too wide for readr's default buffer.
#'
#' @param counts_path Path to the CSV (\code{.csv} or \code{.csv.gz}).
#' @param strip_sample_suffix Logical. If \code{TRUE}, drop a trailing
#'   \code{.<suffix>} from each sample ID (e.g. \code{SAMPLE.1} -> \code{SAMPLE}),
#'   matching the pipeline's rowname cleanup. Default \code{FALSE}.
#' @return A data frame with samples as rows and positions as columns.
#' @keywords internal
read_counts_csv <- function(counts_path, strip_sample_suffix = FALSE) {
    if (!file.exists(counts_path)) {
        stop("Count matrix not found: ", counts_path)
    }

    # these headers are tens of thousands of columns wide; readr's default
    # connection buffer cannot hold one
    wanted <- 131072 * 30
    old <- Sys.getenv("VROOM_CONNECTION_SIZE", unset = NA)
    if (is.na(old) || suppressWarnings(as.numeric(old)) < wanted || is.na(suppressWarnings(as.numeric(old)))) {
        Sys.setenv(VROOM_CONNECTION_SIZE = format(wanted, scientific = FALSE))
        on.exit(
            if (is.na(old)) Sys.unsetenv("VROOM_CONNECTION_SIZE")
            else Sys.setenv(VROOM_CONNECTION_SIZE = old),
            add = TRUE
        )
    }

    counts <- readr::read_csv(counts_path, show_col_types = FALSE, progress = FALSE)

    if (!"Sample_ID" %in% colnames(counts)) {
        stop("Count matrix must have a 'Sample_ID' column: ", counts_path)
    }

    ids <- counts[["Sample_ID"]]
    counts <- as.data.frame(counts[, setdiff(colnames(counts), "Sample_ID"), drop = FALSE])
    if (strip_sample_suffix) {
        ids <- sub("\\.[^.]*$", "", ids)
    }
    rownames(counts) <- ids
    counts
}


#' Run Torino end to end for one gene, from a count CSV to an RDS
#'
#' Reads a per-gene count matrix, preprocesses it with \code{prepare_counts},
#' factorizes it with \code{torino_factorization}, and writes the result to an
#' RDS in the layout the downstream Torino tooling expects.
#'
#' The saved list has three elements:
#' \describe{
#'   \item{gene}{The gene name.}
#'   \item{torino_K<K>}{The \code{torino_factorization} result, named for the
#'     requested \code{K} (e.g. \code{torino_K10}), or \code{NULL} if the
#'     factorization failed.}
#'   \item{coords_all}{All position column names from the original matrix,
#'     before any filtering, column merging or padding.}
#' }
#'
#' @param gene_name Character. Gene name, stored as \code{$gene} in the output.
#' @param K Integer. Number of factors. Also names the \code{torino_K<K>}
#'   element of the output.
#' @param counts_path Path to the count CSV (\code{.csv} or \code{.csv.gz}) with
#'   a \code{Sample_ID} column; remaining columns are genomic positions.
#' @param out_path Path of the RDS to write. Parent directories are created if
#'   needed. Pass \code{NULL} to skip writing and only return the result.
#' @param min_counts Integer. Minimum total read count per sample
#'   (default: 30).
#' @param n Number of consecutive position columns to merge, or the string
#'   \code{"auto"} to pick one from the matrix width via
#'   \code{choose_merge_n}. Default 1 (no merging).
#' @param seed Integer random seed for the factorization (default: 1).
#' @param strip_sample_suffix Logical. If \code{TRUE}, drop a trailing
#'   \code{.<suffix>} from each sample ID before filtering (default
#'   \code{FALSE}).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#' @return Invisibly, the list that was saved.
#' @examples
#' \dontrun{
#' run_torino_gene(
#'     gene_name   = "SRSF3",
#'     K           = 10,
#'     counts_path = "example_data/Human_SRSF3.csv.gz",
#'     out_path    = "results/rds/Human_SRSF3.rds"
#' )
#' }
#' @export
run_torino_gene <- function(gene_name, K, counts_path, out_path,
                            min_counts = 30, n = 1, seed = 1,
                            strip_sample_suffix = FALSE, verbose = TRUE) {

    K <- as.integer(K)
    if (is.na(K) || K < 1) stop("K must be a positive integer.")

    say <- function(...) if (verbose) message(sprintf(...))

    say("Loading %s ...", counts_path)
    counts <- read_counts_csv(counts_path, strip_sample_suffix = strip_sample_suffix)

    coords_all   <- colnames(counts)
    nbases_total <- ncol(counts)
    say("Loaded %d samples x %d bases.", nrow(counts), nbases_total)

    if (identical(as.character(n), "auto")) {
        n <- choose_merge_n(nbases_total)
        say("Column merging: n = %d (auto, from %d bases).", n, nbases_total)
    } else {
        n <- as.integer(n)
        if (is.na(n) || n < 1) stop("n must be a positive integer or \"auto\".")
    }

    say("Preparing counts (min_counts = %d, n = %d) ...", min_counts, n)
    prepared <- prepare_counts(counts, min_counts = min_counts, n = n)
    say("Prepared matrix: %d samples x %d positions (%d samples dropped).",
        nrow(prepared), ncol(prepared), nrow(counts) - nrow(prepared))

    say("Running torino_factorization (K = %d, seed = %d) ...", K, seed)
    fit <- torino_factorization(prepared, K = K, seed = seed)

    if (is.null(fit)) {
        # torino_factorization swallows errors and returns NULL
        warning("torino_factorization returned NULL (factorization failed) for ",
                gene_name)
    } else {
        say("Factorization complete. EF_smooth %d x %d | EL %d x %d | final ELBO %.1f",
            nrow(fit$EF_smooth), ncol(fit$EF_smooth),
            nrow(fit$EL), ncol(fit$EL), utils::tail(fit$elbo, 1))
    }

    result <- list(gene = gene_name)
    result[[paste0("torino_K", K)]] <- fit
    result$coords_all <- coords_all

    if (!is.null(out_path)) {
        dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
        saveRDS(result, file = out_path)
        say("Saved %s (%.1f MB).", out_path, file.info(out_path)$size / 1024^2)
    }

    invisible(result)
}
