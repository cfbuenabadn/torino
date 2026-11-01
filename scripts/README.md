# `scripts/` — deprecated, kept for the pipeline

**Nothing in this directory is part of the installable `torino` R package.**
`.Rbuildignore` excludes `^scripts$`, so none of it ships or installs. The
`scripts/` directory you see in an installed `torino` is `inst/scripts/`,
which is a different thing.

As of torino 0.2.0 this directory is **deprecated**. It is kept — not deleted —
because `torino_paper` and other downstream code load these files by
hardcoded absolute path.