# =============================================================================
# 02_DEG_DESeq2.R
# DESeq2 Wald-test differential expression for both bulk contrasts.
# Aβ42 contrast is fully specified; BKO reuses the same code block
# (template-level, as it comes from a prior in-house project).
# =============================================================================

source(here("02_bulkRNAseq", "00_config.R"))

run_deg <- function(comp) {
  cat("\n=== DE:", comp$name, "===\n")
  dds <- readRDS(file.path(result_dir, paste0("bulk_", comp$name, "_dds.rds")))

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "treatment", "control"),
                 alpha = deseq_padj_cutoff)

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_symbol") %>%
    arrange(padj)

  deg <- res_df %>%
    filter(!is.na(padj), padj < deseq_padj_cutoff,
           abs(log2FoldChange) >= deseq_lfc_cutoff)
  deg_up <- deg %>% filter(log2FoldChange > 0)
  deg_dn <- deg %>% filter(log2FoldChange < 0)

  out_dir <- file.path(result_dir, paste0("bulk_", comp$name, "_DEG"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  write.csv(res_df, file.path(out_dir, paste0(comp$name, "_all_results.csv")),
            row.names = FALSE)
  write.csv(deg,    file.path(out_dir, paste0(comp$name, "_DEG.csv")),
            row.names = FALSE)
  write.csv(deg_up, file.path(out_dir, paste0(comp$name, "_UP.csv")),
            row.names = FALSE)
  write.csv(deg_dn, file.path(out_dir, paste0(comp$name, "_DOWN.csv")),
            row.names = FALSE)

  cat(sprintf("  Genes tested: %d | DEGs: %d (up %d, down %d)\n",
              nrow(res_df), nrow(deg), nrow(deg_up), nrow(deg_dn)))
  invisible(res_df)
}

for (comp in contrasts) {
  run_deg(comp)
}
