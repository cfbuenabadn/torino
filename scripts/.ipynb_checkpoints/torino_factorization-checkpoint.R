.libPaths(c('/project/yangili1/cfbuenabadn/R/x86_64-pc-linux-gnu-library/4.1',
            '/software/R-4.1.0-el8-x86_64/lib64/R/library')) 

library(tidyverse)
library(ebpmf.alpha)

torino_factorization <- function(counts, K, seed=0) {

    torino_out <- tryCatch(
        {
        lib_size <- rowSums(counts) %>% as.integer()
        set.seed(seed)
        fit = ebpmf.alpha::ebpmf_identity(as.matrix(counts),K, lib_size=lib_size)

        factors <- paste0('factor_', 1:K)
        
        rownames(fit$EF_smooth) <- colnames(counts)
        colnames(fit$EF_smooth) <- factors

        rownames(fit$EL) <- rownames(counts)
        colnames(fit$EL) <- factors
            
    
        torino_out <- list(EF = fit$EF,
                           EF_smooth = fit$EF_smooth,
                           EL = fit$EL,
                           elbo = fit$elbo,
                           coords = colnames(counts),
                           samples = rownames(counts)
                          )
        },
        error = function(e) {
            torino_out <- NULL
            }
        )

    return(torino_out)
    
    }


group_matrix <- function(df, n = 3){
  n_groups <- ncol(df) %/% n
  
  # Split dataframe into groups of n columns and sum them
  output <- lapply(seq_len(n_groups), function(i) {
    group_start_col <- n * (i - 1) + 1
    group_end_col <- min(n * i, ncol(df))
    
    group_sum <- rowSums(df[, group_start_col:group_end_col], na.rm = TRUE)
    
    # Create a new data frame with one column and set its name
    group_df <- data.frame(group_sum)
    colnames(group_df) <- colnames(df)[group_start_col]
    
    return(group_df)
  })
  
  # Combine the output into a single dataframe
  output_df <- do.call(cbind, output)
  
  return(output_df)
}

prepare_counts <- function(counts, min_counts, n=1){
    counts <- counts[rowSums(counts) >= min_counts, ]
    
    if (n > 1) {
        counts <- group_matrix(counts, n)
    }
    
    kept_samples <- (counts %>% dim())[1]
    coords <- colnames(counts)
    
    matrix_cols <- c('pois1', coords, 'pois2')
    nsamples <- dim(counts)[1]
    
    set.seed(1)
    counts <- cbind(rpois(nsamples, 0.2), counts, rpois(nsamples, 0.2)) %>% as.matrix()
    return(counts)
    
}

run_torino <- function(counts, K = 10, min_counts = 100){
    counts <- prepare_counts(counts, min_counts)
    torino_factorization <- run_torino(counts, K=K)
    return(torino_fit)
}
