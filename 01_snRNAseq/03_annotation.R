# =============================================================================
# 03_annotation.R
# Cell-type annotation of both Original and DecontX objects with the same
# marker panel and the same cluster -> cell-type map.
# Astrocytes = clusters 1, 7, 8 in both objects.
# =============================================================================

source(here("01_snRNAseq", "00_config.R"))

original <- readRDS(file.path(result_dir, "01_Original_Integrated.rds"))
decontx  <- readRDS(file.path(result_dir, "02_DecontX_Integrated.rds"))

# --- Marker panel ------------------------------------------------------------
marker_panel <- list(
  Astrocyte      = c("Slc1a2", "Slc1a3", "Aqp4", "Sox9", "Gja1", "Aldh1l1"),
  Microglia      = c("Itgam", "Aif1", "Cx3cr1", "Tyrobp", "Trem2", "Hexb", "P2ry12"),
  Neuron         = c("Snap25", "Syt1", "Rbfox3", "Map2"),
  Oligodendrocyte= c("Mog", "Mbp", "Plp1", "Mag"),
  OPC            = c("Pdgfra", "Cspg4"),
  Endothelial    = c("Cldn5", "Pecam1", "Flt1"),
  Pericyte       = c("Rgs5", "Pdgfrb", "Acta2")
)

# --- Cluster -> cell-type map (19 clusters -> 7 types) ---------------------
# Astrocytes: clusters 1, 7, 8
cluster_celltype_map <- c(
  "0" = "Microglia",
  "1" = "Astrocyte",
  "2" = "Neuron",
  "3" = "Oligodendrocyte",
  "4" = "Microglia",
  "5" = "Neuron",
  "6" = "OPC",
  "7" = "Astrocyte",
  "8" = "Astrocyte",
  "9" = "Oligodendrocyte",
  "10" = "Endothelial",
  "11" = "Microglia",
  "12" = "Neuron",
  "13" = "Astrocyte",
  "14" = "Oligodendrocyte",
  "15" = "Pericyte",
  "16" = "Neuron",
  "17" = "Microglia",
  "18" = "Astrocyte"
)

annotate_object <- function(obj, map, name) {
  DefaultAssay(obj) <- "RNA"
  obj$cell_type <- map[as.character(obj$seurat_clusters)]
  obj$cell_type[is.na(obj$cell_type)] <- "Unassigned"

  # Marker dot plot
  markers_use <- unlist(marker_panel)
  markers_use <- intersect(markers_use, rownames(obj))
  p_dot <- DotPlot(obj, features = markers_use, group.by = "cell_type") +
    RotatedAxis() + ggtitle(paste0("Marker profile (", name, ")"))
  ggsave(file.path(figure_dir, paste0("03_", name, "_marker_dotplot.pdf")),
         p_dot, width = 14, height = 6)

  # UMAP coloured by cell type
  p_umap <- DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE)
  ggsave(file.path(figure_dir, paste0("03_", name, "_celltype_umap.pdf")),
         p_umap, width = 8, height = 6)

  saveRDS(obj, file.path(result_dir, paste0("03_", name, "_Annotated.rds")))
  obj
}

original_annotated <- annotate_object(original, cluster_celltype_map, "Original")
decontx_annotated  <- annotate_object(decontx,  cluster_celltype_map, "DecontX")
