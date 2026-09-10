# =============================================================================
# 07_within_snRNA_comparison.R
# Compare the three snRNA astrocyte versions (V1, V2, V3):
#   - Immune-term shrinkage with progressive correction (111 -> 26 -> 4)
#   - Conserved non-immune programs shared by all three versions
#     (metabolism, transport, adhesion/cytoskeleton — all same direction)
# =============================================================================

source(here("01_snRNAseq", "00_config.R"))

all_terms <- read.csv(
  file.path(result_dir, "06_GSEA_AllVersions_WithFunctionalGroups.csv"),
  stringsAsFactors = FALSE
)

# --- 1. Immune-term shrinkage ------------------------------------------------
immune_terms <- all_terms %>%
  filter(Functional_Group == "Immune response") %>%
  group_by(Version) %>%
  summarise(n_immune_terms = n_distinct(ID), .groups = "drop")

immune_terms$Version <- factor(
  immune_terms$Version,
  levels = c("V1_Original", "V2_DecontX", "V3_DecontX_Micro")
)
immune_terms <- immune_terms[order(immune_terms$Version), ]

write.csv(
  immune_terms,
  file.path(result_dir, "07_ImmuneTerms_Shrinkage_ThreeVersions.csv"),
  row.names = FALSE
)
cat("Immune GO-BP terms (padj<0.05, |NES|>=1.5):\n")
print(immune_terms)

# --- 2. Conserved non-immune programs (shared by all 3 versions) ------------
conserved_groups <- c(
  "Metabolism / biosynthesis",
  "Transport / trafficking",
  "Adhesion / cytoskeleton"
)

shared_nonimmune <- all_terms %>%
  filter(Functional_Group %in% conserved_groups) %>%
  group_by(Functional_Group, ID, Description) %>%
  summarise(
    n_versions   = n_distinct(Version),
    n_positive   = sum(NES > 0),
    n_negative   = sum(NES < 0),
    all_same_dir = (n_positive == n_versions | n_negative == n_versions),
    .groups = "drop"
  ) %>%
  filter(n_versions == 3)

# Tally by functional group
conserved_tally <- shared_nonimmune %>%
  group_by(Functional_Group) %>%
  summarise(
    n_shared_terms     = n_distinct(ID),
    n_same_direction   = sum(all_same_dir),
    .groups = "drop"
  )

write.csv(
  shared_nonimmune,
  file.path(result_dir, "07_ConservedNonImmune_SharedByAllThree.csv"),
  row.names = FALSE
)
write.csv(
  conserved_tally,
  file.path(result_dir, "07_ConservedNonImmune_Tally.csv"),
  row.names = FALSE
)
cat("\nConserved non-immune programs (shared by V1, V2, V3):\n")
print(conserved_tally)

# --- 3. Version-unique terms --------------------------------------------------
version_unique <- all_terms %>%
  group_by(ID, Description, Functional_Group) %>%
  summarise(
    n_versions = n_distinct(Version),
    versions   = paste(sort(unique(Version)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_versions == 1)

write.csv(
  version_unique,
  file.path(result_dir, "07_VersionUnique_Terms.csv"),
  row.names = FALSE
)
