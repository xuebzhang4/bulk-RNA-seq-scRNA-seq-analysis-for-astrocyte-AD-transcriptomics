# =============================================================================
# 05_DEG_GSEA_three_versions.R
# Differential expression and gene-set enrichment analysis (GSEA) for each of
# the three astrocyte versions (V1, V2, V3), using an identical engine:
#
#   DE   : Seurat FindMarkers (Wilcoxon, FAD vs WT, logfc.threshold = 0,
#          min.pct = 0.05, RNA assay)
#   Rank : genes ranked by avg_log2FC
#   GSEA : clusterProfiler::gseGO (GO-BP/MF/CC, org.Mm.eg.db,
#          minGSSize = 15, maxGSSize = 150, BH adjustment, by = "fgsea")
#
# Reported threshold (relaxed): padj < 0.05 & |NES| >= 1.5.
# A stricter cut (0.01 / 2.0) is also exported for reference.
# =============================================================================

source(here("01_snRNAseq", "00_config.R"))

v1 <- readRDS(file.path(result_dir, "04_V1_Original_Astrocytes.rds"))
v2 <- readRDS(file.path(result_dir, "04_V2_DecontX_Astrocytes.rds"))
v3 <- readRDS(file.path(result_dir, "04_V3_DecontX_Micro_Astrocytes.rds"))

astro_objects <- list(
  V1_Original       = v1,
  V2_DecontX        = v2,
  V3_DecontX_Micro  = v3
)

# --- Core DE + GSEA pipeline -------------------------------------------------
run_deg_gsea <- function(obj, analysis_name) {
  cat("\n=== DE + GSEA:", analysis_name, "===\n")
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- "orig.ident"
  stopifnot(all(c("FAD", "WT") %in% levels(Idents(obj))))

  # Differential expression (FAD vs WT)
  deg <- FindMarkers(
    obj,
    ident.1        = "FAD",
    ident.2        = "WT",
    logfc.threshold = 0,
    min.pct        = 0.05,
    test.use       = "wilcox",
    verbose        = FALSE
  )
  deg$gene_symbol <- rownames(deg)
  fc_col <- if ("avg_log2FC" %in% colnames(deg)) "avg_log2FC" else "avg_logFC"
  deg$rank_metric <- deg[[fc_col]]

  write.csv(
    deg,
    file.path(result_dir, paste0("05_", analysis_name, "_DEG_FAD_vs_WT.csv")),
    row.names = FALSE
  )

  # Ranked gene list for GSEA (by avg_log2FC)
  gsea_input <- deg %>%
    filter(!is.na(gene_symbol), gene_symbol != "", is.finite(rank_metric)) %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
    ungroup()
  gene_list <- sort(
    setNames(gsea_input$rank_metric, gsea_input$gene_symbol),
    decreasing = TRUE
  )

  # GSEA for each GO ontology
  run_gsego <- function(ontology) {
    gsea_obj <- gseGO(
      geneList     = gene_list,
      OrgDb        = org.Mm.eg.db,
      keyType      = "SYMBOL",
      ont          = ontology,
      minGSSize    = gsea_min_gs_size,
      maxGSSize    = gsea_max_gs_size,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      verbose      = FALSE,
      seed         = TRUE,
      by           = "fgsea"
    )
    all_res <- as.data.frame(gsea_obj)
    if (nrow(all_res) == 0) return(list(all = data.frame(), relaxed = data.frame()))

    all_res$Analysis  <- analysis_name
    all_res$Ontology  <- ontology

    relaxed <- all_res %>%
      filter(p.adjust < gsea_padj_relaxed, abs(NES) >= gsea_nes_relaxed) %>%
      arrange(desc(abs(NES)))

    write.csv(
      all_res,
      file.path(result_dir, paste0("05_", analysis_name, "_GO_", ontology, "_ALL.csv")),
      row.names = FALSE
    )
    write.csv(
      relaxed,
      file.path(result_dir, paste0("05_", analysis_name, "_GO_", ontology, "_RELAXED.csv")),
      row.names = FALSE
    )
    list(all = all_res, relaxed = relaxed)
  }

  go_bp <- run_gsego("BP")
  go_mf <- run_gsego("MF")
  go_cc <- run_gsego("CC")

  # Summary
  data.frame(
    Analysis       = analysis_name,
    Cells_Total    = ncol(obj),
    Cells_WT       = sum(obj$orig.ident == "WT"),
    Cells_FAD      = sum(obj$orig.ident == "FAD"),
    Genes_Ranked   = length(gene_list),
    GO_BP_Relaxed  = nrow(go_bp$relaxed),
    GO_MF_Relaxed  = nrow(go_mf$relaxed),
    GO_CC_Relaxed  = nrow(go_cc$relaxed),
    stringsAsFactors = FALSE
  )
}

# --- Run for all three versions ----------------------------------------------
all_summaries <- lapply(names(astro_objects), function(nm) {
  run_deg_gsea(astro_objects[[nm]], nm)
})
summary_df <- bind_rows(all_summaries)

write.csv(
  summary_df,
  file.path(result_dir, "05_GSEA_Summary_Three_Versions.csv"),
  row.names = FALSE
)
print(summary_df)
