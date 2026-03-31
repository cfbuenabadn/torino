#' Run Torino factorization on a count matrix
#'
#' Decomposes a (samples x positions) read-count matrix into K factors using
#' Empirical Bayes Poisson Matrix Factorization (\code{ebpmf_identity}).
#'
#' @param counts Numeric matrix (rows = samples, columns = genomic positions).
#'   Typically produced by \code{prepare_counts}.
#' @param K Integer. Number of factors to estimate.
#' @param seed Integer random seed (default: 0).
#' @return A named list with components:
#'   \describe{
#'     \item{EF}{Raw factor matrix (positions x K).}
#'     \item{EF_smooth}{Smoothed factor matrix (positions x K).}
#'     \item{EL}{Loading matrix (samples x K).}
#'     \item{elbo}{Evidence lower bound trace.}
#'     \item{coords}{Column names of the input count matrix.}
#'     \item{samples}{Row names of the input count matrix.}
#'   }
#'   Returns \code{NULL} if the factorization fails.
#' @export
torino_factorization <- function(counts, K, seed = 0) {
    torino_out <- tryCatch(
        {
            lib_size <- rowSums(counts) %>% as.integer()
            set.seed(seed)
            fit <- ebpmf_identity(as.matrix(counts), K, lib_size = lib_size)

            factors <- paste0("factor_", seq_len(K))

            rownames(fit$EF_smooth) <- colnames(counts)
            colnames(fit$EF_smooth) <- factors

            rownames(fit$EL) <- rownames(counts)
            colnames(fit$EL) <- factors

            list(
                EF       = fit$EF,
                EF_smooth = fit$EF_smooth,
                EL       = fit$EL,
                elbo     = fit$elbo,
                coords   = colnames(counts),
                samples  = rownames(counts)
            )
        },
        error = function(e) NULL
    )
    return(torino_out)
}


#' Merge consecutive columns by summing them in groups
#'
#' Reduces the width of a count matrix by summing every \code{n} consecutive
#' columns into one. Any trailing columns that do not fill a complete group are
#' dropped. Useful for large genomic windows where RAM is a bottleneck.
#'
#' @param df A data frame or matrix with named columns.
#' @param n Integer. Number of consecutive columns to merge (sum) into one.
#' @return A data frame with \code{ncol(df) \%/\% n} columns. Each column is
#'   named after the first column in its group.
#' @export
group_matrix <- function(df, n = 3) {
    n_groups <- ncol(df) %/% n
    output <- lapply(seq_len(n_groups), function(i) {
        group_start_col <- n * (i - 1) + 1
        group_end_col   <- min(n * i, ncol(df))
        group_sum <- rowSums(df[, group_start_col:group_end_col, drop = FALSE],
                             na.rm = TRUE)
        group_df <- data.frame(group_sum)
        colnames(group_df) <- colnames(df)[group_start_col]
        return(group_df)
    })
    do.call(cbind, output)
}


#' Prepare a count matrix for Torino factorization
#'
#' Filters samples by minimum total counts, optionally merges consecutive
#' columns to reduce RAM usage, and flanks the matrix with Poisson noise
#' boundary columns (lambda = 0.2) to stabilise edge effects during
#' factorization.
#'
#' @param counts Numeric matrix or data frame (rows = samples,
#'   columns = genomic positions).
#' @param min_counts Integer. Minimum total read count per sample; samples
#'   below this threshold are removed.
#' @param n Integer. Number of consecutive columns to merge via
#'   \code{group_matrix} (default: 1, no merging). Typical values for large
#'   genes: 3 (>15 k bases), 9 (>60 k bases), 30 (>1 M bases).
#' @return A numeric matrix ready for \code{torino_factorization}, with two
#'   Poisson-noise boundary columns prepended and appended.
#' @export
prepare_counts <- function(counts, min_counts, n = 1) {
    counts <- counts[rowSums(counts) >= min_counts, , drop = FALSE]

    if (n > 1) {
        counts <- group_matrix(counts, n)
    }

    nsamples <- nrow(counts)
    set.seed(1)
    counts <- cbind(rpois(nsamples, 0.2), counts, rpois(nsamples, 0.2))
    as.matrix(counts)
}


#' End-to-end Torino wrapper
#'
#' Convenience function that calls \code{prepare_counts} followed by
#' \code{torino_factorization}.
#'
#' @param counts Numeric matrix or data frame (rows = samples,
#'   columns = genomic positions).
#' @param K Integer. Number of factors (default: 10).
#' @param min_counts Integer. Minimum total read count per sample
#'   (default: 100).
#' @param n Integer. Number of consecutive columns to merge for RAM reduction
#'   (default: 1, no merging).
#' @param seed Integer random seed (default: 0).
#' @return Result of \code{torino_factorization}, or \code{NULL} on error.
#' @export
run_torino <- function(counts, K = 10, min_counts = 100, n = 1, seed = 0) {
    counts <- prepare_counts(counts, min_counts, n = n)
    torino_factorization(counts, K = K, seed = seed)
}
