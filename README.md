# Astrocyte responses to Aβ pathology — snRNA-seq × bulk RNA-seq analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.4.3-blue.svg)]()
[![Seurat](https://img.shields.io/badge/Seurat-4.3.0-8856a7.svg)]()
[![DESeq2](https://img.shields.io/badge/DESeq2-1.46.0-4daf4a.svg)]()
[![DOI](https://img.shields.io/badge/Zenodo-archive%20pending-lightgrey.svg)]()

Reproducible computational pipeline for comparing astrocyte transcriptional responses
to amyloid-β (Aβ) pathology across **single-nucleus RNA-seq** of the 5xFAD cortex and
two **bulk RNA-seq** experiments on purified primary astrocytes (direct Aβ42 exposure and
astrocytic *Bmal1* knockout).

The central methodological design is to resolve one snRNA-seq dataset into **three
nested astrocyte versions** — Original, DecontX-corrected, and DecontX with residual
microglial-signal nuclei removed — and then place every dataset under **one shared GSEA
/ functional-group framework** so that terms and leading-edge genes can be matched by GO
ID and compared by direction.

**Analysis workflow:** [`docs/workflow/analysis_workflow_flowchart.html`](docs/workflow/analysis_workflow_flowchart.html)
(also available as an editable PowerPoint:
[`docs/workflow/analysis_workflow_flowchart.pptx`](docs/workflow/analysis_workflow_flowchart.pptx)).

---

## Datasets

| | Modality / system | Contrast | Source |
|---|---|---|---|
| **snRNA** | single-nucleus, mouse cortex | WT vs 5xFAD, saline arms only | GEO [**GSE227157**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227157) (Fatmi *et al.*, *Aging* 2024;16:3137–3159, PMID 38385967) |
| **bulk Aβ42** | bulk, purified primary astrocytes | Aβ42 48 h vs untreated, n=3 vs 3 | in-house — *GEO accession to be added* |
| **bulk BKO** | bulk, purified astrocytes | *Bmal1*-KO vs WT, n=3 vs 4 | GEO [**GSE325658**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325658) (prior in-house project) |

Exact sample lists are in [`data/metadata/`](data/metadata). Large raw/count files are
**not** tracked here — see [`data/README.md`](data/README.md) for download instructions.

## Repository structure

```
.
├── 01_snRNAseq/        # MAIN LINE (detailed): QC → DecontX → annotation →
│                        #   three astrocyte versions → DE+GSEA → functional grouping
├── 02_bulkRNAseq/      # DESeq2 for Aβ42 and BKO contrasts (concise; BKO template-level)
├── 03_cross_dataset/   # Cross-modality rules: GO-ID match, NES direction, shared leading edge
├── 04_figures/         # Reusable plotting templates + unified theme/palette
├── data/metadata/      # Small sample sheets (raw data downloaded separately)
├── docs/workflow/      # Workflow diagram (html / pptx)
├── results/            # Outputs (only small final tables committed)
├── README.md
├── LICENSE
├── CITATION.cff
└── .gitignore
```

Each numbered folder has its own README with a file-by-file plan and the intended level
of detail.

## Analysis at a glance

### 1. snRNA-seq main line

1. **QC and CCA integration** — `CreateSeuratObject` (min.cells=3, min.features=200);
   QC filters `200 < nFeature < 6000`, `percent.mt < 20`; LogNormalize; 2000 HVG (vst);
   CCA integration (dims 1:30); PCA (30); FindNeighbors (dims 1:10); clustree sweep;
   final resolution 0.3; UMAP (dims 1:10). → **Original object**.
2. **DecontX correction** — raw counts re-read per sample; `celda::decontX` with
   `z` = Original clusters; decontXcounts extracted; identical QC → CCA → res 0.3 →
   UMAP pipeline. → **DecontX object**.
3. **Annotation** — same marker panel and same cluster→cell-type map for both objects
   (19 clusters → 7 types); **astrocytes = clusters 1, 7, 8**.
4. **Three nested astrocyte versions**:
   - **V1 Original** — clusters 1/7/8 from the Original object (raw counts)
   - **V2 DecontX** — clusters 1/7/8 from the DecontX object (corrected counts)
   - **V3 DecontX-Micro** — built **only from V2**: AddModuleScore of 13 microglial
     genes, remove nuclei with score ≥ 90th percentile. V1 is never modified.
5. **DE + GSEA** — identical engine on each version: `FindMarkers` (Wilcoxon, FAD vs WT,
   logfc.threshold=0, min.pct=0.05); rank by avg_log2FC; `clusterProfiler::gseGO`
   (GO-BP/MF/CC, org.Mm.eg.db, minGSSize=15, maxGSSize=150, BH, by="fgsea").
   Reported threshold: `padj < 0.05 & |NES| ≥ 1.5` (353 / 156 / 87 GO-BP terms).
6. **Functional grouping** — keyword-regex dictionary → 10 groups; align versions by GO
   ID and NES direction.
7. **Within-snRNA comparison** — immune-term shrinkage 111 → 26 → 4; conserved
   non-immune programs (metabolism, transport, adhesion/cytoskeleton — all same direction).

### 2. Bulk RNA-seq

- **DESeq2** Wald DE for Aβ42 vs untreated (n=3 vs 3) and Bmal1-KO vs WT (n=3 vs 4).
- Genes ranked by the **DESeq2 Wald statistic** (deliberately different from snRNA's
  avg_log2FC); same GSEA settings and functional grouping (911 terms for Aβ42; 308 for
  BKO).

### 3. Cross-dataset comparison

- Terms matched by **GO ID**; NES sign compared; shared leading-edge genes = intersection
  of `core_enrichment` strings.
- Per-modality log2FC used for **direction only** — never compared as effect size across
  platforms.
- Combined heatmaps capped at **±1.5**.

## Reproducibility

**Software** (R 4.4.3):

| Package | Version | Purpose |
|---|---|---|
| Seurat | 4.3.0 | snRNA-seq QC, integration, clustering, DE |
| celda | 1.22.1 | DecontX ambient-RNA correction |
| DESeq2 | 1.46.0 | Bulk RNA-seq differential expression |
| clusterProfiler | 4.14.6 | GSEA (gseGO, fgsea) |
| org.Mm.eg.db | 3.20.0 | Mouse gene annotation |
| SummarizedExperiment | 1.36.0 | SCE for DecontX |
| ggplot2 / pheatmap / patchwork | — | Plotting |

All paths are relative and centralised in each folder's `00_config.R` — no personal-drive
paths. GSEA is reported at `padj < 0.05 & |NES| ≥ 1.5` (GO biological process;
BH-adjusted; fgsea implementation; minGSSize 15 / maxGSSize 150).

## Quick start

```r
# 1. Download data into data/raw/ following data/README.md
# 2. From the repository root, run scripts in order:

# snRNA-seq main line
source("01_snRNAseq/00_config.R")          # then 01 -> 07

# Bulk RNA-seq
source("02_bulkRNAseq/00_config.R")        # then 01 -> 04

# Cross-dataset comparison
source("03_cross_dataset/00_shared_utils.R")  # then 01 -> 04
```
