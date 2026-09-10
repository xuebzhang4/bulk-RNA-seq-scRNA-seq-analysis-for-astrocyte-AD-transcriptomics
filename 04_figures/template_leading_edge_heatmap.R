# =============================================================================
# template_leading_edge_heatmap.R
# Reusable pheatmap template for leading-edge gene x dataset log2FC heatmaps.
#
# Expected input:
#   mat: numeric matrix (genes x datasets) of log2FC values, already capped.
#        Use make_capped_lfc_matrix() from 03_cross_dataset/00_shared_utils.R
#        to build and cap the matrix.
# =============================================================================

source(here("04_figures", "theme_palette.R"))
source(here("03_cross_dataset", "00_shared_utils.R"))

# --- Leading-edge heatmap (capped log2FC) ------------------------------------
plot_le_heatmap <- function(mat, title = "", cap = 1.5,
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             file = NULL, width = 5, height = NULL) {
  if (is.null(height)) height <- max(4, 0.22 * nrow(mat) + 1.5)

  pheatmap(
    mat,
    color         = heatmap_palette(100),
    breaks        = seq(-cap, cap, length.out = 101),
    cluster_rows  = cluster_rows,
    cluster_cols  = cluster_cols,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    angle_col     = "45",
    main          = title,
    border_color  = "white",
    na_col        = "grey90",
    fontsize_row  = 10,
    fontsize_col  = 11,
    filename      = file,
    width         = width,
    height        = height
  )
  invisible(mat)
}

# --- Build matrix from a list of log2FC vectors ------------------------------
# lfc_list: named list of named vectors (gene -> log2FC), one per dataset.
# genes: character vector of genes to include (in order).
build_le_matrix <- function(lfc_list, genes, cap = 1.5) {
  lfc_df <- bind_log2fc(lfc_list, gene_universe = genes)
  make_capped_lfc_matrix(lfc_df, cap = cap, gene_order = genes)
}

# --- Example usage (commented out) -------------------------------------------
# lfc_list <- list(
#   V3     = setNames(deg_v3$avg_log2FC, deg_v3$gene_symbol),
#   AB42   = setNames(deg_ab$log2FoldChange, deg_ab$gene_symbol),
#   BKO    = setNames(deg_bko$log2FoldChange, deg_bko$gene_symbol)
# )
# genes <- c("Hmgcr", "Sqle", "Dhcr7", "Hmgcs1", "Fdft1")
# mat <- build_le_matrix(lfc_list, genes, cap = 1.5)
# plot_le_heatmap(mat, title = "Sterol biosynthesis leading-edge genes",
#                 file = "sterol_le_heatmap.pdf")
