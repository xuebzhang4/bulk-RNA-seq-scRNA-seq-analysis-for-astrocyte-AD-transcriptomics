# =============================================================================
# template_umap_dotplot.R
# Reusable templates for UMAP and marker dot/feature plots.
#
# Expected input:
#   seurat_obj: a Seurat object with a "umap" reduction and cell-type labels
#               in seurat_obj$cell_type (or another metadata column).
# =============================================================================

source(here("04_figures", "theme_palette.R"))

# --- UMAP coloured by cell type ----------------------------------------------
plot_umap_celltype <- function(seurat_obj, group_by = "cell_type",
                                label = TRUE, title = "") {
  DimPlot(seurat_obj, reduction = "umap", group.by = group_by,
          label = label, repel = TRUE) +
    scale_color_manual(values = palette_celltypes) +
    ggtitle(title) +
    theme_publication() +
    theme(legend.title = element_blank())
}

# --- UMAP split by condition --------------------------------------------------
plot_umap_split <- function(seurat_obj, split_by = "orig.ident",
                             group_by = "cell_type", title = "") {
  DimPlot(seurat_obj, reduction = "umap", group.by = group_by,
          split.by = split_by) +
    scale_color_manual(values = palette_celltypes) +
    ggtitle(title) +
    theme_publication()
}

# --- Marker dot plot ----------------------------------------------------------
plot_marker_dotplot <- function(seurat_obj, marker_list, group_by = "cell_type",
                                 title = "") {
  markers_use <- unlist(marker_list)
  markers_use <- intersect(markers_use, rownames(seurat_obj))
  DotPlot(seurat_obj, features = markers_use, group.by = group_by) +
    RotatedAxis() +
    ggtitle(title) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

# --- Feature plot (expression on UMAP) ---------------------------------------
plot_feature <- function(seurat_obj, gene, title = NULL, cols = c("grey90", "firebrick3")) {
  if (is.null(title)) title <- gene
  FeaturePlot(seurat_obj, features = gene, reduction = "umap",
              cols = cols, order = TRUE, pt.size = 0.5) +
    ggtitle(title) +
    theme_publication()
}

# --- Example usage (commented out) -------------------------------------------
# p1 <- plot_umap_celltype(seurat_obj, title = "Cell types")
# ggsave("umap_celltypes.pdf", p1, width = 8, height = 6)
#
# markers <- list(Astrocyte = c("Slc1a2","Aqp4"), Microglia = c("Itgam","Trem2"))
# p2 <- plot_marker_dotplot(seurat_obj, markers, title = "Marker expression")
# ggsave("marker_dotplot.pdf", p2, width = 12, height = 5)
