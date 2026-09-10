# =============================================================================
# template_GSEA_nes_bar.R
# Reusable templates for GSEA NES bar charts and functional-group stacked bars.
#
# Expected input:
#   gsea_df: data.frame with columns ID, Description, NES, p.adjust,
#            Functional_Group (optional), Dataset (optional).
# =============================================================================

source(here("04_figures", "theme_palette.R"))

# --- Horizontal NES bar chart (top terms) ------------------------------------
plot_nes_bar <- function(gsea_df, n_top = 20, title = "",
                          padj_cutoff = 0.05, nes_cutoff = 1.5) {
  df <- gsea_df %>%
    filter(p.adjust < padj_cutoff, abs(NES) >= nes_cutoff) %>%
    arrange(desc(abs(NES))) %>%
    head(n_top) %>%
    mutate(Description = factor(Description, levels = rev(Description)))

  ggplot(df, aes(x = NES, y = Description, fill = NES > 0)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("FALSE" = palette_direction["Negative"],
                                  "TRUE"  = palette_direction["Positive"]),
                      labels = c("FALSE" = "Down", "TRUE" = "Up")) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(title = title, x = "Normalised Enrichment Score (NES)",
         y = NULL, fill = "Direction") +
    theme_publication()
}

# --- Functional-group tally bar chart -----------------------------------------
plot_group_tally <- function(gsea_df, group_col = "Functional_Group",
                               dataset_col = "Dataset", title = "") {
  df <- gsea_df %>%
    filter(!is.na(.data[[group_col]])) %>%
    group_by(.data[[group_col]], .data[[dataset_col]]) %>%
    summarise(n_terms = n_distinct(ID), .groups = "drop")

  ggplot(df, aes(x = .data[[group_col]], y = n_terms,
                 fill = .data[[dataset_col]])) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(title = title, x = NULL, y = "Number of GO terms",
         fill = "Dataset") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# --- NES comparison across datasets (dot plot) -------------------------------
plot_nes_dotplot <- function(gsea_df, term_ids, group_col = "Description",
                               dataset_col = "Dataset", title = "") {
  df <- gsea_df %>%
    filter(ID %in% term_ids) %>%
    mutate(Description = factor(Description, levels = unique(Description)))

  ggplot(df, aes(x = .data[[dataset_col]], y = Description,
                  size = -log10(p.adjust), colour = NES)) +
    geom_point() +
    scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3",
                           midpoint = 0) +
    labs(title = title, x = NULL, y = NULL,
         size = "-log10(padj)", colour = "NES") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# --- Example usage (commented out) -------------------------------------------
# p <- plot_nes_bar(gsea_df, n_top = 25, title = "Top GO-BP terms (V3)")
# ggsave("gsea_nes_bar.pdf", p, width = 8, height = 8)
