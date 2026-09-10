# =============================================================================
# 03_BKO_vs_AB42_vs_V3.R
# Three-way comparison: Bmal1-KO (bulk) x Aβ42 (bulk) x V3 (snRNA).
#
#   - Aβ42 ∩ BKO overlap
#   - BKO ∩ V3 but NOT Aβ42 (BMAL1-associated, in-vivo-only programs)
#   - Direction concordance (same-direction vs opposite-direction)
# =============================================================================

source(here("03_cross_dataset", "00_shared_utils.R"))

result_dir <- here("results")

load_bp <- function(prefix) {
  read.csv(file.path(result_dir, paste0(prefix, "_GO_BP_RELAXED.csv")),
           stringsAsFactors = FALSE)
}

gsea_three <- list(
  Bulk_BKO          = load_bp("bulk_BKO_vs_WT"),
  Bulk_AB42         = load_bp("bulk_AB42_vs_untreated"),
  V3_DecontX_Micro  = load_bp("05_V3_DecontX_Micro")
)

aligned <- align_by_goid(gsea_three, padj_cutoff = 0.05, nes_cutoff = 1.5)
wide <- aligned$wide

# --- Aβ42 ∩ BKO --------------------------------------------------------------
ab_bko <- wide %>%
  filter(Significant_Bulk_AB42, Significant_Bulk_BKO) %>%
  mutate(Direction_Concordant =
           (NES_Bulk_AB42 > 0 & NES_Bulk_BKO > 0) |
           (NES_Bulk_AB42 < 0 & NES_Bulk_BKO < 0))

write.csv(ab_bko,
          file.path(result_dir, "03_AB42_intersect_BKO.csv"),
          row.names = FALSE)
cat("Aβ42 ∩ BKO:", nrow(ab_bko), "terms;",
    sum(ab_bko$Direction_Concordant, na.rm = TRUE), "same-direction\n")

# --- BKO ∩ V3 but NOT Aβ42 ---------------------------------------------------
bko_v3_not_ab <- wide %>%
  filter(Significant_Bulk_BKO, Significant_V3_DecontX_Micro,
         !Significant_Bulk_AB42) %>%
  mutate(Direction_Concordant =
           (NES_Bulk_BKO > 0 & NES_V3_DecontX_Micro > 0) |
           (NES_Bulk_BKO < 0 & NES_V3_DecontX_Micro < 0))

write.csv(bko_v3_not_ab,
          file.path(result_dir, "03_BKO_intersect_V3_not_AB42.csv"),
          row.names = FALSE)
cat("BKO ∩ V3 (not Aβ42):", nrow(bko_v3_not_ab), "terms;",
    sum(bko_v3_not_ab$Direction_Concordant, na.rm = TRUE), "same-direction;",
    sum(!bko_v3_not_ab$Direction_Concordant, na.rm = TRUE), "opposite-direction\n")

# --- Shared leading-edge genes for BKO ∩ V3 not Aβ42 ------------------------
bko_v3_genes <- shared_le_table(
  list(Bulk_BKO = gsea_three$Bulk_BKO, V3 = gsea_three$V3_DecontX_Micro),
  padj_cutoff = 0.05, nes_cutoff = 1.5
)
if (!is.null(bko_v3_genes) && nrow(bko_v3_genes) > 0) {
  write.csv(bko_v3_genes,
            file.path(result_dir, "03_BKO_V3_SharedLeadingEdge.csv"),
            row.names = FALSE)
}

# --- Summary ------------------------------------------------------------------
summary_df <- data.frame(
  Comparison = c("Aβ42 ∩ BKO", "BKO ∩ V3 (not Aβ42)"),
  N_Terms    = c(nrow(ab_bko), nrow(bko_v3_not_ab)),
  Same_Direction = c(sum(ab_bko$Direction_Concordant, na.rm = TRUE),
                     sum(bko_v3_not_ab$Direction_Concordant, na.rm = TRUE)),
  Opposite_Direction = c(sum(!ab_bko$Direction_Concordant, na.rm = TRUE),
                         sum(!bko_v3_not_ab$Direction_Concordant, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(summary_df,
          file.path(result_dir, "03_ThreeWay_Comparison_Summary.csv"),
          row.names = FALSE)
print(summary_df)
