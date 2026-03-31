# Torino
Reference-free analysis of mRNA isoform variation in large-scale RNA-seq datasets

## Installation

Torino depends on several packages from GitHub. Install them first, then install torino itself:

```r
install.packages("remotes")

# GitHub dependencies not on CRAN
remotes::install_github("stephenslab/mixsqp")
remotes::install_github("stephenslab/smashr")
remotes::install_github("stephenslab/ebnm")
remotes::install_github("stephenslab/fastTopics")
remotes::install_github("willwerscheid/flashier")
remotes::install_github("DongyueXie/smashrgen")
remotes::install_github("DongyueXie/ebpm")
remotes::install_github("DongyueXie/vebpm")

# Install torino
remotes::install_github("yangili1/cfbuenabadn/torino")
```

CRAN dependencies (`Matrix`, `magrittr`, `wavethresh`, `matrixStats`, `Rfast`, `ashr`) are installed automatically.

> **Note:** The package was developed under R 4.1.0. It has not been tested on newer R versions. Known risk: `Rfast::Pmin` was removed in Rfast 2.x; if you encounter errors in the VGA solver, replace calls to `Pmin()` with base R `pmin()`.

## Usage

### Quick start (R)

```r
library(torino)

# counts: matrix with rows = samples, columns = genomic positions
fit <- run_torino(counts, K = 10, min_counts = 100)

# fit$EF_smooth  — smoothed factor matrix (positions x K)
# fit$EL         — loading matrix (samples x K)
# fit$elbo       — ELBO trace
```

`run_torino()` is an end-to-end wrapper. For more control, call the steps separately:

```r
# 1. Filter samples and optionally merge columns to reduce RAM
#    n=3 for >15k bases, n=9 for >60k bases, n=30 for >1M bases
prepped <- prepare_counts(counts, min_counts = 100, n = 3)

# 2. Factorize
fit <- torino_factorization(prepped, K = 10, seed = 1)
```

### Command-line runner

`inst/scripts/run_torino_factorization.R` provides a command-line interface for cluster use:

```bash
Rscript run_torino_factorization.R <input_csv> <min_counts> <K> <out_rds> [merge_n]
```

- `input_csv`: CSV with a `Sample_ID` column (rows = samples, remaining columns = genomic positions)
- `min_counts`: minimum total reads per sample
- `K`: number of factors
- `out_rds`: output RDS file
- `merge_n`: columns to merge per group — integer, or `auto` (default) for automatic selection based on matrix size

Optional: set `TORINO_FILTER_SAMPLES` to a path containing one sample ID per line to restrict the analysis to a subset of samples.

```bash
TORINO_FILTER_SAMPLES=samples.txt Rscript run_torino_factorization.R gene.csv 100 10 out.rds auto
```

The output RDS contains `$gene`, `$torino` (the factorization result), `$K`, and `$coords_all` (position names before any column merging).
