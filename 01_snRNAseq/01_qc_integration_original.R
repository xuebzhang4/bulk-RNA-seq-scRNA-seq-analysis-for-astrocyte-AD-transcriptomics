# =============================================================================
# 01_qc_integration_original.R
# QC, CCA integration, clustering and UMAP for the Original (uncorrected) object.
# Standard pipeline — key decisions: CCA dims 1:30, clustering resolution 0.3
# (selected via clustree), UMAP dims 1:10.
# =============================================================================

source(here("01_snRNAseq", "00_config.R"))

# --- 1. Read 10x counts and build Seurat objects ----------------------------
sample_dirs <- list(
  WT  = file.path(data_dir, "WT"),
  FAD = file.path(data_dir, "FAD")
)

seurat_list <- lapply(names(sample_dirs), function(sam) {
  counts <- Read10X(data.dir = sample_dirs[[sam]])
  sobj <- CreateSeuratObject(
    counts   = counts,
    project  = sam,
    min.cells    = seurat_min_cells,
    min.features = seurat_min_features
  )
  mt_pattern <- if (any(grepl("^mt-", rownames(sobj)))) "^mt-" else "^MT-"
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = mt_pattern)
  sobj
})
names(seurat_list) <- names(sample_dirs)

# --- 2. QC filtering ---------------------------------------------------------
seurat_list <- lapply(seurat_list, function(x) {
  subset(x, subset = nFeature_RNA > qc_min_features &
                    nFeature_RNA < qc_max_features &
                    percent.mt   < qc_max_percent_mt)
})

# --- 3. Normalisation and variable features ----------------------------------
seurat_list <- lapply(seurat_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n_hvg)
  x
})

# --- 4. CCA integration (sequential to avoid Windows deadlocks) -------------
plan("sequential")
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = n_hvg)
anchors  <- FindIntegrationAnchors(
  object.list    = seurat_list,
  anchor.features = features,
  dims           = cca_dims,
  reduction      = "cca"
)
seurat_original <- IntegrateData(anchorset = anchors, dims = cca_dims)
plan("multisession", workers = 4)

# --- 5. Scale, PCA, neighbours, clustering ----------------------------------
DefaultAssay(seurat_original) <- "integrated"
seurat_original <- ScaleData(seurat_original)
seurat_original <- RunPCA(seurat_original, npcs = 30, approx = FALSE)
seurat_original <- FindNeighbors(seurat_original, dims = pca_dims)

# Resolution sweep for clustree evaluation
seurat_original <- FindClusters(
  seurat_original,
  resolution = c(0.1, 0.2, 0.3, 0.4, 0.5)
)
p_clustree <- clustree(seurat_original, prefix = "integrated_snn_res.")
ggsave(file.path(figure_dir, "01_clustree.pdf"), p_clustree, width = 8, height = 10)

# Final resolution = 0.3
seurat_original <- FindClusters(seurat_original, resolution = cluster_resolution)
Idents(seurat_original) <- "integrated_snn_res.0.3"
seurat_original <- RunUMAP(seurat_original, reduction = "pca", dims = pca_dims)

# --- 6. Save -----------------------------------------------------------------
saveRDS(
  seurat_original,
  file.path(result_dir, "01_Original_Integrated.rds")
)

p_umap <- DimPlot(seurat_original, reduction = "umap", label = TRUE)
ggsave(file.path(figure_dir, "01_Original_UMAP.pdf"), p_umap, width = 7, height = 6)
