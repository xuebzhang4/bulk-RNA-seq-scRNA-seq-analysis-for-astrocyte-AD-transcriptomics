# 01_snRNAseq — single-nucleus RNA-seq main line

This is the **most detailed** part of the repository: the ambient-RNA correction
strategy and the three nested astrocyte versions are the methodological core of the
paper, so scripts `02`–`05` are complete, runnable and heavily commented. Standard
QC/integration (`01`) is kept concise on purpose.

Run order matches the numbering. All scripts read paths from `00_packages_and_config.R`
(relative paths — no hard-coded `D:\...` drives).

| Script | What it does | Detail level |
|---|---|---|
| `00_packages_and_config.R` | Package loading, version pins, relative paths, and every shared parameter in one place (QC cut-offs, CCA/PCA dims, resolution, GSEA thresholds, the 13 microglia genes). | Support — full |
| `01_qc_integration_original.R` | Build the **Original** object: `CreateSeuratObject(min.cells=3, min.features=200)`; QC `200 < nFeature < 6000`, `percent.mt < 20`; `LogNormalizeData`; 2000 HVG (`vst`); CCA integration (`dims 1:30`); `ScaleData`; PCA (30); `FindNeighbors(dims 1:10)`; clustree sweep; final `resolution = 0.3`; UMAP (`dims 1:10`). | **Condensed** — standard pipeline, short comments; only the decisions that matter (resolution choice via clustree, dims) are spelled out. |
| `02_decontX_corrected.R` | Build the parallel **DecontX** object: re-read raw counts per sample, run `celda::decontX` **with `z` = the Original clusters from `01`**, extract `decontXcounts`, rebuild a Seurat object and repeat the *identical* QC → CCA → res 0.3 → UMAP path. | **Full** — the `z`-from-cluster-1 design and per-sample execution are made explicit because they are easy to get wrong. |
| `03_celltype_annotation.R` | Annotate **both** objects with the same marker panel and the same cluster→cell-type map (19 clusters → 7 types). Astrocytes = clusters **1 / 7 / 8** in both objects. | Medium — marker list and cluster map are shown in full; the plotting is minimal. |
| `04_astrocyte_three_versions.R` | ★ Derive the three nested astrocyte objects: **V1 Original** = clusters 1/7/8 from `01` (raw RNA); **V2 DecontX** = clusters 1/7/8 from `02` (corrected counts); **V3 DecontX-Micro** is built **only from V2** — `AddModuleScore` of 13 microglial genes (P2ry12, Tmem119, Hexb, Cx3cr1, Aif1, Itgam, C1qa/b/c, Trem2, Tyrobp, Ptprc, Lgals3), then remove nuclei with score ≥ the 90th percentile (top 10 %). An explicit comment states V1 is never modified. | **Most detailed** — this cascade is the paper's key control and must be reproducible line by line. |
| `05_DEG_GSEA_three_versions.R` | ★ Run the identical engine on each version: `FindMarkers(5xFAD vs WT, test.use="wilcox", logfc.threshold=0, min.pct=0.05, assay="RNA")`; rank genes by **avg_log2FC**; `clusterProfiler::gseGO` (BP/MF/CC, `org.Mm.eg.db`, `minGSSize=15, maxGSSize=150`, BH, `by="fgsea"`). Report the relaxed set `padj<0.05 & |NES|≥1.5` (353 / 156 / 87); the strict 0.01 / 2.0 cut is kept as a commented option. | **Most detailed.** |
| `06_functional_grouping.R` | Collapse GO terms into the 10 functional groups (and finer subgroups) with a keyword-regex dictionary; align the three versions **by GO ID**, compare NES direction, tag leading-edge genes that carry microglia/DAM markers. The manual curation/renaming step is clearly separated from the automatic one. | Medium. |
| `07_within_snRNA_comparison.R` | Version-to-version comparison: immune-term shrinkage 111 → 26 → 4; conserved non-immune programs shared by all three (metabolism 21, transport 13, adhesion/cytoskeleton 12, all same-direction); subgroup and leading-edge tables. | Medium. |

## Inputs / outputs

- **Input:** `data/raw/gse227157/` 10x outputs + `data/metadata/snRNA_samples.csv`.
- **Output (git-ignored, regenerated):** integrated Original / DecontX objects and the
  three astrocyte `.rds`, per-version DE and GSEA tables.
- **Small final tables** that are curated/renamed for the paper may be copied to
  `../results/tables/` and committed.
