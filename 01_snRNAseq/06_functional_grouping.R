# =============================================================================
# 06_functional_grouping.R
# Collapse GO terms into functional groups using a keyword-regex dictionary,
# then align the three versions by GO ID and NES direction.
# Automatic grouping is followed by a manual curation step (clearly marked).
# =============================================================================

source(here("01_snRNAseq", "00_config.R"))

# --- Load GSEA results for all three versions --------------------------------
load_gsea_relaxed <- function(version) {
  read.csv(
    file.path(result_dir, paste0("05_", version, "_GO_BP_RELAXED.csv")),
    stringsAsFactors = FALSE
  )
}

gsea_list <- list(
  V1_Original      = load_gsea_relaxed("V1_Original"),
  V2_DecontX       = load_gsea_relaxed("V2_DecontX"),
  V3_DecontX_Micro = load_gsea_relaxed("V3_DecontX_Micro")
)

# --- Functional-group keyword dictionary (10 groups) ------------------------
# Terms are assigned by regex matching against Description + ID.
# The first matching group wins; unmatched terms go to "Other".
group_patterns <- list(
  "Immune response" = c("immune", "inflamm", "antigen", "MHC", "complement",
                         "cytokine", "interferon", "leukocyte", "lymphocyte",
                         "defense response", "tolerogen", "immunoglobulin"),
  "Metabolism / biosynthesis" = c("metabolic", "biosynthesis", "catabolic",
                         "oxidative phosphorylation", "respiratory chain",
                         "ATP synthesis", "glycolysis", "tricarboxylic",
                         "lipid metabolic", "steroid", "sterol", "cholesterol",
                         "fatty acid", "amino acid", "nucleotide", "NAD",
                         "mitochondrial", "electron transport"),
  "Transport / trafficking" = c("transport", "trafficking", "vesicle",
                         "endocytosis", "exocytosis", "secretion", "import",
                         "export", "glutamate", "amino acid transport",
                         "ion transport", "ER to Golgi", "intraciliary",
                         "extracellular fluid"),
  "Adhesion / cytoskeleton" = c("cytoskeleton", "microtubule", "actin",
                         "intermediate filament", "cell adhesion", "focal adhesion",
                         "axoneme", "cilium", "flagellum", "cell junction",
                         "extracellular matrix", "basement membrane"),
  "Signaling" = c("signal transduction", "signaling pathway", "MAPK", "Wnt",
                         "Notch", "TGF-beta", "BMP", "Hippo", "PI3K", "AKT",
                         "mTOR", "Ras", "GPCR", "second messenger"),
  "Cell cycle / proliferation" = c("cell cycle", "mitotic", "DNA replication",
                         "chromosome segregation", "spindle", "proliferation",
                         "apoptosis", "cell death", "autophagy"),
  "Transcription / chromatin" = c("transcription", "chromatin", "histone",
                         "RNA splicing", "RNA processing", "ribosome",
                         "translation", "epigenetic", "methylation", "acetylation"),
  "Development / differentiation" = c("development", "differentiation",
                         "morphogenesis", "patterning", "neurogenesis",
                         "gliogenesis", "axonogenesis", "synapse"),
  "Stress response" = c("stress response", "unfolded protein", "heat shock",
                         "oxidative stress", "reactive oxygen", "DNA damage",
                         "xenobiotic", "detoxification"),
  "Cellular organisation" = c("cellular component", "organelle", "cellular component assembly",
                         "cellular component organisation", "cellular component biogenesis")
)

assign_functional_group <- function(description) {
  desc_lower <- tolower(description)
  for (grp in names(group_patterns)) {
    if (any(sapply(group_patterns[[grp]], function(p) grepl(tolower(p), desc_lower)))) {
      return(grp)
    }
  }
  "Other"
}

# --- Apply automatic grouping to each version --------------------------------
for (nm in names(gsea_list)) {
  gsea_list[[nm]]$Functional_Group <-
    sapply(gsea_list[[nm]]$Description, assign_functional_group)
  gsea_list[[nm]]$Version <- nm
}

# --- MANUAL CURATION (override automatic assignments where needed) -----------
# Add specific GO ID -> group overrides here. Example:
# manual_overrides <- c("GO:0030317" = "Other")
manual_overrides <- c()

for (nm in names(gsea_list)) {
  idx <- match(gsea_list[[nm]]$ID, names(manual_overrides))
  gsea_list[[nm]]$Functional_Group[!is.na(idx)] <- manual_overrides[idx[!is.na(idx)]]
}

# --- Align three versions by GO ID -------------------------------------------
all_terms <- bind_rows(gsea_list)
all_terms$Direction <- ifelse(all_terms$NES > 0, "Positive", "Negative")

# Wide format: one row per GO term, columns for each version's NES
nes_wide <- all_terms %>%
  select(ID, Description, Functional_Group, Version, NES, Direction) %>%
  pivot_wider(names_from = Version, values_from = c(NES, Direction), names_sep = "_")

write.csv(
  all_terms,
  file.path(result_dir, "06_GSEA_AllVersions_WithFunctionalGroups.csv"),
  row.names = FALSE
)
write.csv(
  nes_wide,
  file.path(result_dir, "06_GSEA_AlignedByGOID_NES_Direction.csv"),
  row.names = FALSE
)

# --- Functional-group tally ---------------------------------------------------
group_tally <- all_terms %>%
  group_by(Functional_Group, Version) %>%
  summarise(n_terms = n_distinct(ID), .groups = "drop") %>%
  pivot_wider(names_from = Version, values_from = n_terms, values_fill = 0)

write.csv(
  group_tally,
  file.path(result_dir, "06_FunctionalGroup_Tally_ThreeVersions.csv"),
  row.names = FALSE
)
print(group_tally)
