# =============================================================================
# 04_culture_purity_check.R
# Brief check of astrocyte culture purity in the Aβ42 bulk data:
# expression of canonical astrocyte, microglia, neuron, and oligodendrocyte
# markers across samples, confirming the purified-astrocyte starting material.
# =============================================================================

source(here("02_bulkRNAseq", "00_config.R"))

comp <- contrasts$AB42_vs_untreated
counts <- read.csv(comp$count_file, row.names = 1, check.names = FALSE)

purity_markers <- list(
  Astrocyte      = c("Slc1a2", "Slc1a3", "Aqp4", "Sox9", "Gja1", "Aldh1l1", "Gfap"),
  Microglia      = c("Itgam", "Aif1", "Cx3cr1", "Trem2", "Hexb", "P2ry12"),
  Neuron         = c("Snap25", "Syt1", "Rbfox3", "Map2"),
  Oligodendrocyte= c("Mog", "Mbp", "Plp1")
)

marker_expr <- lapply(names(purity_markers), function(ct) {
  genes <- intersect(purity_markers[[ct]], rownames(counts))
  data.frame(
    CellType = ct,
    Gene     = genes,
    MeanCount = rowMeans(counts[genes, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
})
marker_df <- bind_rows(marker_expr)

write.csv(
  marker_df,
  file.path(result_dir, "bulk_AB42_culture_purity_markers.csv"),
  row.names = FALSE
)

p <- ggplot(marker_df, aes(x = CellType, y = log10(MeanCount + 1), fill = CellType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5) +
  labs(title = "Culture purity marker expression (Aβ42 bulk astrocytes)",
       x = NULL, y = "log10(mean count + 1)") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(figure_dir, "bulk_AB42_culture_purity.pdf"), p, width = 6, height = 5)

cat("Astrocyte markers detected:", sum(marker_df$CellType == "Astrocyte" & marker_df$MeanCount > 0),
    "/", length(purity_markers$Astrocyte), "\n")
cat("Microglia markers detected:", sum(marker_df$CellType == "Microglia" & marker_df$MeanCount > 0),
    "/", length(purity_markers$Microglia), "\n")
