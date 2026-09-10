# =============================================================================
# 00_config.R — shared configuration for 02_bulkRNAseq
# All paths are relative to the repository root. Source this file first.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(here)
})

set.seed(123)

# --- Paths --------------------------------------------------------------------
repo_root  <- here::here()
data_dir   <- file.path(repo_root, "data", "raw")
result_dir <- file.path(repo_root, "results")
figure_dir <- file.path(repo_root, "results", "figures")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# --- DE parameters ------------------------------------------------------------
deseq_padj_cutoff <- 0.05
deseq_lfc_cutoff  <- 0.5
min_count_filter   <- 10   # genes kept if rowSums(counts >= 10) >= n_control_samples

# --- GSEA parameters (identical to snRNA for cross-dataset comparability) ---
gsea_min_gs_size <- 15
gsea_max_gs_size <- 150
gsea_padj_relaxed <- 0.05
gsea_nes_relaxed  <- 1.5

# --- Contrasts ----------------------------------------------------------------
# Only the two astrocyte contrasts used in the manuscript are defined here.
contrasts <- list(
  AB42_vs_untreated = list(
    name            = "AB42_vs_untreated",
    count_file      = file.path(data_dir, "bulk_AB42", "Astrocyte_Only_Clean_Counts.csv"),
    metadata_file   = file.path(repo_root, "data", "metadata", "bulk_AB42_samples.csv"),
    control_samples = c("Astro_1", "Astro_2", "Astro_3"),
    treat_samples   = c("A_Beta_1", "A_Beta_2", "A_Beta_3")
  ),
  BKO_vs_WT = list(
    name            = "BKO_vs_WT",
    count_file      = file.path(data_dir, "bulk_BKO", "BKO_WT_Astrocyte_Counts.csv"),
    metadata_file   = file.path(repo_root, "data", "metadata", "bulk_BKO_samples.csv"),
    control_samples = c("Astro_WT_1", "Astro_WT_2", "Astro_WT_3", "Astro_WT_4"),
    treat_samples   = c("Astro_BKO_1", "Astro_BKO_2", "Astro_BKO_3")
  )
)
