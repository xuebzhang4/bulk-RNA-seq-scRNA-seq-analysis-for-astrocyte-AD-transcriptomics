# =============================================================================
# 01_unified_GO_tally.R
# Tally significant GO-BP terms and functional groups across all four datasets
# (snRNA V1/V2/V3, bulk Aβ42, bulk BKO) under one identical threshold
# (padj < 0.05 & |NES| >= 1.5).
# =============================================================================

source(here("03_cross_dataset", "00_shared_utils.R"))
source(here("01_snRNAseq", "00_config.R"))

result_dir <- here("results")

# --- Load relaxed GO-BP results for all datasets -----------------------------
load_gsea_bp <- function(prefix) {
  read.csv(file.path(result_dir, paste0(prefix, "_GO_BP_RELAXED.csv")),
           stringsAsFactors = FALSE)
}

gsea_all <- list(
  V1_Original       = load_gsea_bp("05_V1_Original"),
  V2_DecontX        = load_gsea_bp("05_V2_DecontX"),
  V3_DecontX_Micro  = load_gsea_bp("05_V3_DecontX_Micro"),
  Bulk_AB42         = load_gsea_bp("bulk_AB42_vs_untreated"),
  Bulk_BKO          = load_gsea_bp("bulk_BKO_vs_WT")
)

# --- Align by GO ID -----------------------------------------------------------
aligned <- align_by_goid(gsea_all, padj_cutoff = 0.05, nes_cutoff = 1.5)

write.csv(
  aligned$wide,
  file.path(result_dir, "03_AllDatasets_AlignedByGOID.csv"),
  row.names = FALSE
)

# --- Tally by dataset ---------------------------------------------------------
dataset_tally <- data.frame(
  Dataset = names(gsea_all),
  GO_BP_Terms = sapply(gsea_all, nrow),
  stringsAsFactors = FALSE
)
write.csv(
  dataset_tally,
  file.path(result_dir, "03_GO_BP_Tally_AllDatasets.csv"),
  row.names = FALSE
)
print(dataset_tally)

# --- Functional-group tally (if Functional_Group column exists) --------------
if ("Functional_Group" %in% colnames(aligned$long)) {
  group_tally <- aligned$long %>%
    filter(Significant) %>%
    group_by(Functional_Group, Dataset) %>%
    summarise(n_terms = n_distinct(ID), .groups = "drop") %>%
    pivot_wider(names_from = Dataset, values_from = n_terms, values_fill = 0)
  write.csv(
    group_tally,
    file.path(result_dir, "03_FunctionalGroup_Tally_AllDatasets.csv"),
    row.names = FALSE
  )
  print(group_tally)
}
