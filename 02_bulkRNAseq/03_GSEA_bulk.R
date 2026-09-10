# =============================================================================
# 03_GSEA_bulk.R
# GSEA for both bulk contrasts.
#
# Key difference from snRNA: genes are ranked by the DESeq2 Wald STATISTIC
# (not avg_log2FC). All other GSEA parameters are identical to the snRNA
# pipeline (minGSSize 15, maxGSSize 150, BH, fgsea, relaxed padj<0.05 &
# |NES|>=1.5) so the two modalities are directly comparable.
# =============================================================================

source(here("02_bulkRNAseq", "00_config.R"))

run_bulk_gsea <- function(comp) {
  cat("\n=== GSEA:", comp$name, "===\n")

  deg_file <- file.path(result_dir, paste0("bulk_", comp$name, "_DEG"),
                         paste0(comp$name, "_all_results.csv"))
  res_df <- read.csv(deg_file, stringsAsFactors = FALSE)

  # Rank by Wald statistic
  gsea_input <- res_df %>%
    filter(!is.na(gene_symbol), gene_symbol != "", is.finite(stat)) %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  gene_list <- sort(
    setNames(gsea_input$stat, gsea_input$gene_symbol),
    decreasing = TRUE
  )

  write.csv(
    data.frame(gene_symbol = names(gene_list), wald_stat = gene_list),
    file.path(result_dir, paste0("bulk_", comp$name, "_ranked_gene_list.csv")),
    row.names = FALSE
  )

  # GSEA for each GO ontology
  run_gsego <- function(ontology) {
    gsea_obj <- gseGO(
      geneList      = gene_list,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = ontology,
      minGSSize     = gsea_min_gs_size,
      maxGSSize     = gsea_max_gs_size,
      pAdjustMethod = "BH",
      pvalueCutoff  = 1,
      verbose       = FALSE,
      seed          = TRUE,
      by            = "fgsea"
    )
    all_res <- as.data.frame(gsea_obj)
    if (nrow(all_res) == 0) return(list(all = data.frame(), relaxed = data.frame()))

    all_res$Analysis <- comp$name
    all_res$Ontology <- ontology

    relaxed <- all_res %>%
      filter(p.adjust < gsea_padj_relaxed, abs(NES) >= gsea_nes_relaxed) %>%
      arrange(desc(abs(NES)))

    write.csv(all_res,
              file.path(result_dir, paste0("bulk_", comp$name, "_GO_", ontology, "_ALL.csv")),
              row.names = FALSE)
    write.csv(relaxed,
              file.path(result_dir, paste0("bulk_", comp$name, "_GO_", ontology, "_RELAXED.csv")),
              row.names = FALSE)
    list(all = all_res, relaxed = relaxed)
  }

  go_bp <- run_gsego("BP")
  go_mf <- run_gsego("MF")
  go_cc <- run_gsego("CC")

  summary_row <- data.frame(
    Analysis      = comp$name,
    Genes_Ranked  = length(gene_list),
    GO_BP_Relaxed = nrow(go_bp$relaxed),
    GO_MF_Relaxed = nrow(go_mf$relaxed),
    GO_CC_Relaxed = nrow(go_cc$relaxed),
    stringsAsFactors = FALSE
  )
  cat(sprintf("  GO-BP relaxed: %d | GO-MF: %d | GO-CC: %d\n",
              nrow(go_bp$relaxed), nrow(go_mf$relaxed), nrow(go_cc$relaxed)))
  summary_row
}

all_summaries <- lapply(contrasts, run_bulk_gsea)
summary_df <- bind_rows(all_summaries)
write.csv(summary_df,
          file.path(result_dir, "bulk_GSEA_Summary_BothContrasts.csv"),
          row.names = FALSE)
print(summary_df)
