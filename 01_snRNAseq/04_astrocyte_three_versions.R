# =============================================================================
# 04_astrocyte_three_versions.R
# Derive three nested astrocyte objects:
#
#   V1  Original        = astrocyte clusters (1,7,8) from the Original object
#                          (raw counts, no ambient-RNA correction)
#   V2  DecontX         = astrocyte clusters (1,7,8) from the DecontX object
#                          (celda::decontX-corrected counts)
#   V3  DecontX_Micro   = V2 with nuclei carrying a high microglia-module score
#                          removed (top 10% by AddModuleScore of 13 microglia
#                          genes). V3 is built ONLY from V2 — V1 is never
#                          modified at this step.
#
# This cascade is the core methodological design of the study.
# =============================================================================

source(here("01_snRNAseq", "00_config.R"))

original <- readRDS(file.path(result_dir, "03_Original_Annotated.rds"))
decontx  <- readRDS(file.path(result_dir, "03_DecontX_Annotated.rds"))

# --- Helper: extract astrocyte clusters and renormalise ---------------------
extract_astrocytes <- function(obj, clusters, version_name) {
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- "seurat_clusters"

  sub <- subset(obj, idents = clusters)
  sub <- NormalizeData(sub, normalization.method = "LogNormalize", verbose = FALSE)
  sub$version <- version_name
  sub
}

# --- V1: Original astrocytes (raw counts) -----------------------------------
v1_original <- extract_astrocytes(
  obj       = original,
  clusters  = astrocyte_clusters,
  version_name = "V1_Original"
)

# --- V2: DecontX astrocytes (corrected counts) ------------------------------
v2_decontx <- extract_astrocytes(
  obj       = decontx,
  clusters  = astrocyte_clusters,
  version_name = "V2_DecontX"
)

# --- V3: DecontX + microglia-score filter (built ONLY from V2) -------------
# NOTE: V1 is never touched here. The filter is applied exclusively to V2.

DefaultAssay(v2_decontx) <- "RNA"
v2_for_filter <- NormalizeData(v2_decontx, normalization.method = "LogNormalize", verbose = FALSE)

# Microglia module score (13 canonical microglia genes)
microglia_genes_use <- intersect(microglia_score_genes, rownames(v2_for_filter))
stopifnot(length(microglia_genes_use) >= 4)

v2_for_filter <- AddModuleScore(
  v2_for_filter,
  features = list(microglia_genes_use),
  name     = "MicrogliaScore"
)
score_col <- grep("^MicrogliaScore", colnames(v2_for_filter@meta.data), value = TRUE)
score_col <- score_col[length(score_col)]
v2_for_filter$MicrogliaScore_Final <- v2_for_filter@meta.data[[score_col]]

# Threshold = 90th percentile (remove top 10% highest microglia-score nuclei)
score_threshold <- quantile(
  v2_for_filter$MicrogliaScore_Final,
  probs = microglia_score_quantile,
  na.rm = TRUE
)
v2_for_filter$High_Microglia_Score <-
  v2_for_filter$MicrogliaScore_Final >= score_threshold

# Diagnostic: microglia-score distribution by condition
p_vln <- VlnPlot(
  v2_for_filter,
  features = "MicrogliaScore_Final",
  group.by = "orig.ident",
  pt.size  = 0.1
) + ggtitle("Microglia module score in DecontX astrocytes (V2)")
ggsave(file.path(figure_dir, "04_V2_microglia_score_violin.pdf"),
       p_vln, width = 6, height = 5)

# Remove high-microglia-score nuclei -> V3
v3_decontx_micro <- subset(
  v2_for_filter,
  subset = High_Microglia_Score == FALSE
)
v3_decontx_micro <- NormalizeData(
  v3_decontx_micro,
  normalization.method = "LogNormalize",
  verbose = FALSE
)
v3_decontx_micro$version <- "V3_DecontX_Micro"

# --- Cell-count summary ------------------------------------------------------
count_summary <- data.frame(
  Version = c("V1_Original", "V2_DecontX", "V3_DecontX_Micro"),
  Total   = c(ncol(v1_original), ncol(v2_decontx), ncol(v3_decontx_micro)),
  WT      = c(sum(v1_original$orig.ident == "WT"),
              sum(v2_decontx$orig.ident == "WT"),
              sum(v3_decontx_micro$orig.ident == "WT")),
  FAD     = c(sum(v1_original$orig.ident == "FAD"),
              sum(v2_decontx$orig.ident == "FAD"),
              sum(v3_decontx_micro$orig.ident == "FAD")),
  stringsAsFactors = FALSE
)
write.csv(
  count_summary,
  file.path(result_dir, "04_astrocyte_three_versions_cell_counts.csv"),
  row.names = FALSE
)
print(count_summary)

# --- Save all three objects --------------------------------------------------
saveRDS(v1_original,       file.path(result_dir, "04_V1_Original_Astrocytes.rds"))
saveRDS(v2_decontx,        file.path(result_dir, "04_V2_DecontX_Astrocytes.rds"))
saveRDS(v3_decontx_micro,  file.path(result_dir, "04_V3_DecontX_Micro_Astrocytes.rds"))
