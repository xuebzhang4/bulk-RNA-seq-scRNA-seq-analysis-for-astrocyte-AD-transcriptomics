# =============================================================================
# 02_AB42_vs_snRNA.R
# Cross-modality comparison: direct Aβ42 (bulk) vs 5xFAD astrocytes (snRNA).
# Focus on two modules:
#   - Antigen presentation / MHC activation
#   - Sterol / steroid biosynthesis suppression
# For each module: cross-dataset NES table and shared leading-edge genes.
# =============================================================================

source(here("03_cross_dataset", "00_shared_utils.R"))

result_dir <- here("results")

# --- Load GO-BP results -------------------------------------------------------
load_bp <- function(prefix) {
  read.csv(file.path(result_dir, paste0(prefix, "_GO_BP_RELAXED.csv")),
           stringsAsFactors = FALSE)
}

gsea_compare <- list(
  V3_DecontX_Micro = load_bp("05_V3_DecontX_Micro"),
  V2_DecontX        = load_bp("05_V2_DecontX"),
  V1_Original       = load_bp("05_V1_Original"),
  Bulk_AB42         = load_bp("bulk_AB42_vs_untreated")
)

# --- Module keyword filters ---------------------------------------------------
antigen_keywords <- c("antigen", "MHC", "major histocompatibility",
                       "antigen processing", "immune response-regulating",
                       "leukocyte", "T cell")
sterol_keywords  <- c("sterol", "steroid", "cholesterol", "isoprenoid",
                       "terpenoid", "lipid biosynthesis")

filter_module <- function(df, keywords) {
  desc <- tolower(paste(df$ID, df$Description))
  keep <- rep(FALSE, nrow(df))
  for (kw in keywords) keep <- keep | grepl(tolower(kw), desc)
  df[keep, ]
}

# --- Antigen/MHC module -------------------------------------------------------
antigen_terms <- lapply(gsea_compare, function(d) filter_module(d, antigen_keywords))
antigen_aligned <- align_by_goid(antigen_terms, padj_cutoff = 0.05, nes_cutoff = 1.5)

write.csv(
  antigen_aligned$wide,
  file.path(result_dir, "03_AB42_vs_snRNA_AntigenMHC_NES.csv"),
  row.names = FALSE
)

antigen_shared_le <- shared_le_table(antigen_terms, padj_cutoff = 0.05, nes_cutoff = 1.5)
if (!is.null(antigen_shared_le) && nrow(antigen_shared_le) > 0) {
  write.csv(antigen_shared_le,
            file.path(result_dir, "03_AB42_vs_snRNA_AntigenMHC_SharedLE.csv"),
            row.names = FALSE)
}

# --- Sterol/steroid module ----------------------------------------------------
sterol_terms <- lapply(gsea_compare, function(d) filter_module(d, sterol_keywords))
sterol_aligned <- align_by_goid(sterol_terms, padj_cutoff = 0.05, nes_cutoff = 1.5)

write.csv(
  sterol_aligned$wide,
  file.path(result_dir, "03_AB42_vs_snRNA_SterolSteroid_NES.csv"),
  row.names = FALSE
)

sterol_shared_le <- shared_le_table(sterol_terms, padj_cutoff = 0.05, nes_cutoff = 1.5)
if (!is.null(sterol_shared_le) && nrow(sterol_shared_le) > 0) {
  write.csv(sterol_shared_le,
            file.path(result_dir, "03_AB42_vs_snRNA_SterolSteroid_SharedLE.csv"),
            row.names = FALSE)
}

cat("Antigen/MHC terms (significant in >=2 datasets):",
    sum(antigen_aligned$wide$N_Significant_Datasets >= 2, na.rm = TRUE), "\n")
cat("Sterol/steroid terms (significant in >=2 datasets):",
    sum(sterol_aligned$wide$N_Significant_Datasets >= 2, na.rm = TRUE), "\n")
