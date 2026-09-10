# =============================================================================
# 00_config.R — shared configuration for 01_snRNAseq
# All paths are relative to the repository root. Source this file first.
# =============================================================================

# --- Packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(celda)
  library(SingleCellExperiment)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(patchwork)
  library(clustree)
  library(dplyr)
  library(tidyr)
  library(future)
  library(here)
})

set.seed(1234)

# --- Paths (relative to repo root) ------------------------------------------
repo_root   <- here::here()
data_dir    <- file.path(repo_root, "data", "raw", "gse227157")
result_dir  <- file.path(repo_root, "results")
figure_dir  <- file.path(repo_root, "results", "figures")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# --- Parallel ----------------------------------------------------------------
plan("multisession", workers = 4)
options(future.globals.maxSize = 8 * 1024^3)

# --- QC & integration parameters ---------------------------------------------
qc_min_features   <- 200
qc_max_features   <- 6000
qc_max_percent_mt <- 20
seurat_min_cells  <- 3
seurat_min_features <- 200
n_hvg             <- 2000
cca_dims          <- 1:30
pca_dims          <- 1:10
cluster_resolution <- 0.3

# --- Astrocyte clusters (both Original and DecontX objects) ----------------
astrocyte_clusters <- c("1", "7", "8")

# --- Microglia-score genes (V3 filtering) -----------------------------------
microglia_score_genes <- c(
  "P2ry12", "Tmem119", "Hexb", "Cx3cr1", "Aif1", "Itgam",
  "C1qa", "Trem2", "Tyrobp", "Ptprc", "C1qb", "C1qc", "Lgals3"
)
microglia_score_quantile <- 0.90

# --- GSEA parameters ---------------------------------------------------------
gsea_min_gs_size <- 15
gsea_max_gs_size <- 150
gsea_padj_relaxed <- 0.05
gsea_nes_relaxed  <- 1.5
# Strict cut (kept for reference, not reported in the manuscript):
gsea_padj_strict <- 0.01
gsea_nes_strict  <- 2.0

# --- Version names -----------------------------------------------------------
version_names <- c(
  V1 = "Original",
  V2 = "DecontX",
  V3 = "DecontX_Micro"
)
