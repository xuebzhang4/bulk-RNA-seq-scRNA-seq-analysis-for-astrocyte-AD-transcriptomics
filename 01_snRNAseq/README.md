# 01_snRNAseq — single-nucleus RNA-seq main line

This is the ambient-RNA correction strategy and the three nested astrocyte versions are the methodological core of the
paper, so scripts `02`–`05` are complete, runnable and heavily commented. Standard
QC/integration (`01`) is kept concise on purpose.

## Inputs / outputs

- **Input:** `data/raw/gse227157/` 10x outputs + `data/metadata/snRNA_samples.csv`.
- **Output (git-ignored, regenerated):** integrated Original / DecontX objects and the
  three astrocyte `.rds`, per-version DE and GSEA tables.
- **Small final tables** that are curated/renamed for the paper may be copied to
  `../results/tables/` and committed.
