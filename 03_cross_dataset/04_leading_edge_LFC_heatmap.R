# =============================================================================
# 04_leading_edge_LFC_heatmap.R
# Generic leading-edge-gene x dataset log2FC heatmap.
#
# Gene-level log2FC is modality-specific (DESeq2 for bulk, FindMarkers
# avg_log2FC for snRNA) and used for DIRECTION ONLY — values are capped at
# +/- 1.5 for visual parity across datasets.
# =============================================================================

source(here("03_cross_dataset", "00_shared_utils.R"))

result_dir <- here("results")
figure_dir <- here("results", "figures")

# --- Load per-dataset DE results and extract log2FC --------------------------
load_sc_deg <- function(version) {
  df <- read.csv(file.path(result_dir, paste0("05_", version, "_DEG_FAD_vs_WT.csv")),
                  stringsAsFactors = FALSE)
  setNames(df$avg_log2FC, df$gene_symbol)
}

load_bulk_deg <- function(comp_name) {
  df <- read.csv(file.path(result_dir, paste0("bulk_", comp_name, "_DEG"),
                           paste0(comp_name, "_all_results.csv")),
                  stringsAsFactors = FALSE)
  setNames(df$log2FoldChange, df$gene_symbol)
}

lfc_list <- list(
  V1_Original      = load_sc_deg("V1_Original"),
  V2_DecontX       = load_sc_deg("V2_DecontX"),
  V3_DecontX_Micro = load_sc_deg("V3_DecontX_Micro"),
  Bulk_AB42        = load_bulk_deg("AB42_vs_untreated"),
  Bulk_BKO         = load_bulk_deg("BKO_vs_WT")
)

# --- Example: leading-edge genes for a specific GO term ----------------------
# Replace GO_TERM_ID with the term of interest, or supply a custom gene list.
GO_TERM_ID <- "GO:0006694"  # example: steroid biosynthetic process

gsea_v3 <- read.csv(file.path(result_dir, "05_V3_DecontX_Micro_GO_BP_RELAXED.csv"),
                     stringsAsFactors = FALSE)
term_row <- gsea_v3[gsea_v3$ID == GO_TERM_ID, ]

if (nrow(term_row) > 0 && !is.na(term_row$core_enrichment) && term_row$core_enrichment != "") {
  leading_genes <- trimws(unlist(strsplit(term_row$core_enrichment, "/")))

  lfc_df <- bind_log2fc(lfc_list, gene_universe = leading_genes)
  mat <- make_capped_lfc_matrix(lfc_df, cap = 1.5, gene_order = leading_genes)

  write.csv(
    data.frame(gene_symbol = rownames(mat), mat, check.names = FALSE),
    file.path(result_dir, paste0("03_LE_", GO_TERM_ID, "_log2FC_matrix.csv")),
    row.names = FALSE
  )

  plot_lfc_heatmap(
    mat,
    title = paste0(term_row$Description[1], "\nLeading-edge log2FC (capped +/-1.5)"),
    cap   = 1.5,
    file  = file.path(figure_dir, paste0("03_LE_", GO_TERM_ID, "_heatmap.pdf"))
  )
  cat("Heatmap generated for", GO_TERM_ID, "with", nrow(mat), "genes\n")
} else {
  cat("GO term", GO_TERM_ID, "not found or has no leading-edge genes in V3.\n")
  cat("Supply a custom gene list via gene_universe in bind_log2fc() instead.\n")
}
