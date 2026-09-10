# =============================================================================
# 01_deseq2_qc.R
# Build DESeq2 objects, pre-filter, and run VST-based PCA and sample
# correlation QC for both bulk contrasts (Aβ42 and BKO).
# =============================================================================

source(here("02_bulkRNAseq", "00_config.R"))

build_dds <- function(comp) {
  counts <- read.csv(comp$count_file, row.names = 1, check.names = FALSE)
  meta   <- read.csv(comp$metadata_file, stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample_name

  all_samples <- c(comp$control_samples, comp$treat_samples)
  counts <- counts[, all_samples]
  meta   <- meta[all_samples, ]

  condition <- factor(
    c(rep("control", length(comp$control_samples)),
      rep("treatment", length(comp$treat_samples))),
    levels = c("control", "treatment")
  )
  meta$condition <- condition

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = meta,
    design    = ~ condition
  )
  # Pre-filter: keep genes with >= 10 counts in at least n_control samples
  keep <- rowSums(counts(dds) >= min_count_filter) >= length(comp$control_samples)
  dds <- dds[keep, ]
  dds
}

# --- Build and QC each contrast ----------------------------------------------
for (comp in contrasts) {
  cat("QC:", comp$name, "\n")
  dds <- build_dds(comp)
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE)

  # PCA
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = name), size = 3, max.overlaps = Inf) +
    ggtitle(paste0("PCA: ", comp$name)) +
    theme_minimal()
  ggsave(
    file.path(figure_dir, paste0("bulk_", comp$name, "_PCA.pdf")),
    p_pca, width = 7, height = 6
  )

  # Sample correlation heatmap
  sample_dist <- dist(t(assay(vsd)))
  dist_matrix <- as.matrix(sample_dist)
  pdf(file.path(figure_dir, paste0("bulk_", comp$name, "_correlation.pdf")),
      width = 6, height = 5)
  pheatmap(dist_matrix,
           clustering_distance_rows = sample_dist,
           clustering_distance_cols = sample_dist,
           main = paste0(comp$name, " sample distance"))
  dev.off()

  saveRDS(dds, file.path(result_dir, paste0("bulk_", comp$name, "_dds.rds")))
  cat("  Genes retained:", nrow(dds), "\n")
}
