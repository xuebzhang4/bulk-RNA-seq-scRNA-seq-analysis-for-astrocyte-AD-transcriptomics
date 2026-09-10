# 01_snRNAseq — single-nucleus RNA-seq main line

This is the **most detailed** part of the repository: the ambient-RNA correction
strategy and the three nested astrocyte versions are the methodological core of the
paper, so scripts `02`–`05` are complete, runnable and heavily commented. Standard
QC/integration (`01`) is kept concise on purpose.

Run order matches the numbering. All scripts read paths from `00_packages_and_config.R`
(relative paths — no hard-coded `D:\...` drives).

## Inputs / outputs

- **Input:** `data/raw/gse227157/` 10x outputs + `data/metadata/snRNA_samples.csv`.
- **Output (git-ignored, regenerated):** integrated Original / DecontX objects and the
  three astrocyte `.rds`, per-version DE and GSEA tables.
- **Small final tables** that are curated/renamed for the paper may be copied to
  `../results/tables/` and committed.
