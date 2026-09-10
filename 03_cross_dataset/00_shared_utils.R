# =============================================================================
# 00_shared_utils.R
# Shared functions for cross-dataset comparison. These implement the five
# comparison rules used throughout 03_cross_dataset:
#
#   1. Terms are matched by GO ID (not by term text).
#   2. Direction is compared on NES sign (same-direction vs opposite-direction).
#   3. Shared leading-edge genes = intersection of core_enrichment strings
#      (split on "/") for the same GO term.
#   4. Gene-level log2FC is modality-specific (DESeq2 for bulk, FindMarkers
#      for snRNA) and used for DIRECTION ONLY — never compared as effect size
#      across platforms.
#   5. Combined heatmaps use a fixed colour cap of +/- 1.5.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(pheatmap)
})

# --- Rule 1 + 2: align datasets by GO ID, compare NES direction --------------
# gsea_list: named list of data.frames, each with columns ID, Description, NES,
#            p.adjust, core_enrichment (and optionally Functional_Group)
# padj_cutoff, nes_cutoff: significance threshold applied per dataset
align_by_goid <- function(gsea_list, padj_cutoff = 0.05, nes_cutoff = 1.5) {
  all_terms <- bind_rows(lapply(names(gsea_list), function(nm) {
    d <- gsea_list[[nm]]
    d$Dataset <- nm
    d$Significant <- !is.na(d$p.adjust) & d$p.adjust < padj_cutoff & abs(d$NES) >= nes_cutoff
    d$Direction <- ifelse(d$NES > 0, "Positive", "Negative")
    d
  }))

  # Wide: one row per GO term, NES and direction per dataset
  nes_wide <- all_terms %>%
    select(ID, Description, Dataset, NES, Direction, Significant) %>%
    pivot_wider(names_from = Dataset,
                values_from = c(NES, Direction, Significant),
                names_sep = "_")

  # Count significant datasets and assess direction concordance
  sig_cols <- grep("^Significant_", colnames(nes_wide), value = TRUE)
  nes_cols <- grep("^NES_", colnames(nes_wide), value = TRUE)

  nes_wide$N_Significant_Datasets <- rowSums(nes_wide[, sig_cols], na.rm = TRUE)

  nes_wide$All_Same_Direction <- apply(nes_wide[, nes_cols], 1, function(x) {
    x <- x[!is.na(x)]
    length(x) >= 2 && (all(x > 0) || all(x < 0))
  })

  list(long = all_terms, wide = nes_wide)
}

# --- Rule 3: shared leading-edge genes ---------------------------------------
# core_enrichment strings are split on "/" and intersected per GO term.
shared_leading_edge <- function(gsea_list, go_id) {
  gene_sets <- lapply(gsea_list, function(d) {
    row <- d[d$ID == go_id, ]
    if (nrow(row) == 0 || is.na(row$core_enrichment) || row$core_enrichment == "")
      return(character(0))
    trimws(unlist(strsplit(row$core_enrichment, "/")))
  })
  Reduce(intersect, gene_sets)
}

# Leading-edge genes for all shared terms across datasets
shared_le_table <- function(gsea_list, padj_cutoff = 0.05, nes_cutoff = 1.5) {
  aligned <- align_by_goid(gsea_list, padj_cutoff, nes_cutoff)
  shared_ids <- aligned$wide %>%
    filter(N_Significant_Datasets == length(gsea_list)) %>%
    pull(ID)

  bind_rows(lapply(shared_ids, function(gid) {
    genes <- shared_leading_edge(gsea_list, gid)
    if (length(genes) == 0) return(NULL)
    data.frame(
      ID = gid,
      Description = aligned$wide$Description[aligned$wide$ID == gid][1],
      LeadingEdgeGene = genes,
      stringsAsFactors = FALSE
    )
  }))
}

# --- Rule 4: bind per-dataset log2FC (direction only) ------------------------
# lfc_list: named list of named vectors (gene -> log2FC), one per dataset.
#            Bulk uses DESeq2 log2FoldChange; snRNA uses FindMarkers avg_log2FC.
# Returns a data.frame with one row per gene and one log2FC column per dataset.
bind_log2fc <- function(lfc_list, gene_universe = NULL) {
  if (is.null(gene_universe)) {
    gene_universe <- unique(unlist(lapply(lfc_list, names)))
  }
  df <- data.frame(gene_symbol = gene_universe, stringsAsFactors = FALSE)
  for (nm in names(lfc_list)) {
    df[[paste0("log2FC_", nm)]] <- lfc_list[[nm]][df$gene_symbol]
  }
  df
}

# --- Rule 5: capped log2FC matrix for heatmap --------------------------------
# lfc_df: output of bind_log2fc (gene_symbol + log2FC_* columns)
# cap: fixed colour cap (default 1.5)
# Returns a numeric matrix (genes x datasets), values capped at +/- cap.
make_capped_lfc_matrix <- function(lfc_df, cap = 1.5, gene_order = NULL) {
  lfc_cols <- grep("^log2FC_", colnames(lfc_df), value = TRUE)
  mat <- as.matrix(lfc_df[, lfc_cols, drop = FALSE])
  rownames(mat) <- lfc_df$gene_symbol
  colnames(mat) <- sub("^log2FC_", "", lfc_cols)
  mat[!is.na(mat) & mat >  cap] <-  cap
  mat[!is.na(mat) & mat < -cap] <- -cap
  if (!is.null(gene_order)) mat <- mat[intersect(gene_order, rownames(mat)), , drop = FALSE]
  mat
}

# --- Convenience: pheatmap with fixed diverging palette ----------------------
plot_lfc_heatmap <- function(mat, title = "", cap = 1.5, file = NULL,
                              width = 5, height = NULL) {
  colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  breaks <- seq(-cap, cap, length.out = 101)
  if (is.null(height)) height <- max(4, 0.22 * nrow(mat) + 1.5)
  if (!is.null(file)) {
    pdf(file, width = width, height = height)
    pheatmap(mat, color = colors, breaks = breaks,
             cluster_rows = FALSE, cluster_cols = FALSE,
             show_rownames = TRUE, show_colnames = TRUE,
             angle_col = "45", main = title,
             border_color = "white", na_col = "grey90",
             fontsize_row = 10, fontsize_col = 11)
    dev.off()
  }
  invisible(mat)
}
