# Legacy scripts

These are the original per-step scripts from which the pipeline was developed:
parameter setup (`prior_*`), expression and read simulation (`s01_s02_*`,
`sim_step*`), quantification (`s03_salmon.sh`, `s04_STAR.sh`, `s05_*`, `s06_*`),
QTL file preparation (`s08_*`) and model training (`step4a_*`, `step4b_*`).

They are kept for reference only. They contain paths and queue names from one
cluster, are not maintained, and are not used by anything in `bin/` or `R/`.
Everything they did is now done by `isoqtl-sim` from a configuration file; see the
[README](../README.md) and [docs/GETTING_STARTED.md](../docs/GETTING_STARTED.md).
