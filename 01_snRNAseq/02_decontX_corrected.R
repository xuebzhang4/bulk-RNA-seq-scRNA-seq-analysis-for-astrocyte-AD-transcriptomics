# =============================================================================
# 02_decontX_corrected.R
# Ambient-RNA correction with celda::decontX, followed by the identical
# QC -> CCA -> clustering -> UMAP pipeline used for the Original object.
#
# Key design: decontX is supplied with cluster labels (z) taken from the
# Original integrated object (01), so the contamination model uses the
# already-established cell-type/cluster structure.
# =============================================================================

source(here("01_snRNAseq", "00_config.R"))

# --- 1. Load Original object for cluster labels ------------------------------
original <- readRDS(file.path(result_dir, "01_Original_Integrated.rds"))
cluster_col <- "integrated_snn_res.0.3"
if (!cluster_col %in% colnames(original@meta.data)) {
  cluster_col <- "seurat_clusters"
}

# --- 2. Per-sample decontX correction ----------------------------------------
sample_dirs <- list(
  WT  = file.path(data_dir, "WT"),
  FAD = file.path(data_dir, "FAD")
)

decontx_counts_list <- list()

for (sam in names(sample_dirs)) {
  cat("DecontX correction:", sam, "\n")

  counts <- Read10X(data.dir = sample_dirs[[sam]])
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }

  # Cluster labels from the Original object, matched to this sample's barcodes
  sample_cells <- colnames(original)[original$orig.ident == sam]
  clusters <- as.character(original@meta.data[sample_cells, cluster_col])
  names(clusters) <- sample_cells

  # Align barcodes between Original object and raw 10x matrix
  common_cells <- intersect(colnames(counts), names(clusters))
  counts_use   <- counts[, common_cells, drop = FALSE]
  clusters_use <- as.factor(clusters[common_cells])

  # Run decontX with z = Original clusters
  sce <- SingleCellExperiment(assays = list(counts = counts_use))
  sce <- decontX(x = sce, z = clusters_use)

  corrected <- as(assay(sce, "decontXcounts"), "dgCMatrix")
  decontx_counts_list[[sam]] <- corrected

  cat("  Mean contamination fraction:",
      round(mean(sce$decontX_contamination, na.rm = TRUE), 4), "\n")
}

# --- 3. Rebuild Seurat objects and QC ---------------------------------------
seurat_list_decontx <- lapply(names(decontx_counts_list), function(sam) {
  sobj <- CreateSeuratObject(
    counts   = decontx_counts_list[[sam]],
    project  = sam,
    min.cells    = seurat_min_cells,
    min.features = seurat_min_features
  )
  mt_pattern <- if (any(grepl("^mt-", rownames(sobj)))) "^mt-" else "^MT-"
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = mt_pattern)
  sobj <- RenameCells(sobj, add.cell.id = sam)
  sobj
})
names(seurat_list_decontx) <- names(decontx_counts_list)

# Same QC filters as the Original object
seurat_list_decontx <- lapply(seurat_list_decontx, function(x) {
  subset(x, subset = nFeature_RNA > qc_min_features &
                    nFeature_RNA < qc_max_features &
                    percent.mt   < qc_max_percent_mt)
})

# --- 4. Normalisation, CCA integration (identical to Original) --------------
seurat_list_decontx <- lapply(seurat_list_decontx, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n_hvg)
  x
})

plan("sequential")
features <- SelectIntegrationFeatures(object.list = seurat_list_decontx, nfeatures = n_hvg)
anchors  <- FindIntegrationAnchors(
  object.list     = seurat_list_decontx,
  anchor.features = features,
  dims            = cca_dims,
  reduction       = "cca"
)
seurat_decontx <- IntegrateData(anchorset = anchors, dims = cca_dims)
plan("multisession", workers = 4)

# --- 5. Scale, PCA, clustering (resolution 0.3), UMAP ----------------------
DefaultAssay(seurat_decontx) <- "integrated"
seurat_decontx <- ScaleData(seurat_decontx)
seurat_decontx <- RunPCA(seurat_decontx, npcs = 30, approx = FALSE)
seurat_decontx <- FindNeighbors(seurat_decontx, dims = pca_dims)
seurat_decontx <- FindClusters(seurat_decontx, resolution = cluster_resolution)
Idents(seurat_decontx) <- "integrated_snn_res.0.3"
seurat_decontx <- RunUMAP(seurat_decontx, reduction = "pca", dims = pca_dims)

# --- 6. Save -----------------------------------------------------------------
saveRDS(
  seurat_decontx,
  file.path(result_dir, "02_DecontX_Integrated.rds")
)

p_umap <- DimPlot(seurat_decontx, reduction = "umap", label = TRUE)
ggsave(file.path(figure_dir, "02_DecontX_UMAP.pdf"), p_umap, width = 7, height = 6)
