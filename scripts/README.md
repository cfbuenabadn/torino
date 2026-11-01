# `scripts/` — deprecated, kept for the pipeline

**Nothing in this directory is part of the installable `torino` R package.**
`.Rbuildignore` excludes `^scripts$`, so none of it ships or installs. The
`scripts/` directory you see in an installed `torino` is `inst/scripts/`,
which is a different thing.

As of torino 0.2.0 this directory is **deprecated**. It is kept — not deleted —
because `torino_paper` and other downstream code load these files by
**hardcoded absolute path**, e.g.

```python
sys.path.append('/project2/yangili1/cfbuenabadn/torino/scripts')
from torino_isoforms import *
```

```r
source("/project2/yangili1/cfbuenabadn/torino/scripts/torino_factorization.R")
```

Deleting or moving a file here breaks those callers silently. Treat the
directory as a frozen compatibility surface: fix bugs here only when the
pipeline needs it, and put new work in the canonical location instead.

## Where the current versions live

| file | status | canonical version |
|---|---|---|
| `torino_isoforms.py` | deprecated, still imported | `torino_utils/scripts/torino_isoforms.py` |
| `torino_plots.py` | deprecated, still imported | `torino_utils/scripts/torino_plots.py` |
| `torino_factorization.R` | superseded by the package | `torino::prepare_counts`, `torino::torino_factorization`, `torino::run_torino_gene` |
| `plot_factors.R` | **sole copy in the tree** — not yet migrated | — |
| `run_ebpmf_filtered.R` | dead | — |
| `torino_isoforms_old.py` | dead | — |

`run_ebpmf_filtered.R` and `torino_isoforms_old.py` have no known callers.
`run_ebpmf_filtered.R` is gtex-stm era: it imports `ebpmf.alpha` (the
pre-rename package) and sources
`/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts/process_factor_functions.R`.

`plot_factors.R` (`plot_factors_wrapper`) has no copy anywhere else, so it is
deprecated only in the sense that it sits in a deprecated directory — it is
still the only implementation. It needs a home before this directory can be
retired.

## Known divergence: hardcoded `torino_K10`

The copy of `torino_isoforms.py` here reads the factorization out of an RDS with
a **hardcoded** key:

```python
rds['rds']['torino_K10']['EF_smooth']     # torino/scripts  (this file)
rds[f'torino_K{str(K)}']['EF_smooth']     # torino_utils/scripts  (current)
```

Everything in `torino_paper` imports *this* copy, so RDS files written at
`K != 10` — which `torino::run_torino_gene` makes easy to produce — will fail
here while reading correctly through `torino_utils`. Until the two are
reconciled, assume K = 10 for any consumer that loads from this directory.

## `torino_factorization.R` specifically

Its `torino_factorization`, `group_matrix` and `prepare_counts` are now
exported by the package. The two differ only in cosmetics and dead code —
identical formals, identical arithmetic, verified equal on a test matrix — so
sourcing this file shadows the package functions with equivalent behavior. New
code should `library(torino)` instead. Note that the dead
`matrix_cols <- c('pois1', coords, 'pois2')` line here is never applied, which
is why the two Poisson pad columns are named `rpois(nsamples, 0.2)` in saved
RDS files; the package reproduces that for compatibility.
