# 02_bulkRNAseq — bulk RNA-seq (DESeq2)

Clean and reproducible, but **concise**: only the two astrocyte contrasts used in the
paper are kept — Aβ42 vs untreated and Bmal1-KO vs WT. (The old H2O2 and neuron
contrasts are deliberately removed.) The BKO contrast comes from a prior in-house
project and is kept at **template level** — enough to reproduce the numbers, without
the exploratory code from that project.

The same GSEA settings and functional grouping as `01_snRNAseq` are used so the two
modalities are directly comparable. The one deliberate difference is the GSEA ranking
statistic (see `03`).

| Script | What it does | Detail level |
|---|---|---|
| `00_packages_and_config.R` | Packages, version pins, relative paths, contrast definitions for both experiments. | Support — full |
| `01_deseq2_qc.R` | Read count matrix + sample sheet; `DESeqDataSetFromMatrix(~condition)`; pre-filter genes (≥10 counts in the required number of samples); `vst`; PCA and sample-correlation QC. | Concise — QC is standard, only the pre-filter rule is highlighted. |
| `02_DEG_DESeq2.R` | Wald-test DE for the two contrasts: **Aβ42 vs untreated** (n=3 vs 3) and **BKO vs WT** (n=3 vs 4). Output per-gene log2FC, Wald statistic, adjusted p. LFC shrinkage (`apeglm`) provided as a commented option. | Full for Aβ42; BKO reuses the same code block with a different metadata file (template level). |
| `03_GSEA_bulk.R` | Rank genes by the **DESeq2 Wald statistic** (not avg_log2FC — this is the cross-modality difference) and run `gseGO` with the identical settings as snRNA (`minGSSize=15, maxGSSize=150`, BH, `by="fgsea"`, `padj<0.05 & |NES|≥1.5`) and the same 10-group dictionary. Produces 911 terms (Aβ42) and 308 (BKO). | Medium. |
| `04_culture_purity_check.R` | Short check of astrocyte culture purity markers on the Aβ42 bulk data (confirms the purified-astrocyte starting material). | Brief. |

## Inputs / outputs

- **Input:** `data/raw/bulk_AB42/`, `data/raw/bulk_BKO/` and the sample sheets in
  `../data/metadata/`.
- **Output:** DE and GSEA tables per contrast (regenerated; large objects git-ignored).
