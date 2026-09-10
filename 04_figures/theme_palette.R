# =============================================================================
# theme_palette.R
# Unified ggplot2 theme and colour palette for all figures.
# Source this file in every plotting script.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# --- Publication theme --------------------------------------------------------
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", linewidth = 0.5),
      axis.text        = element_text(colour = "black"),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      legend.title     = element_text(face = "bold"),
      legend.position  = "right",
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text       = element_text(face = "bold")
    )
}

# --- Palette: three snRNA versions -------------------------------------------
palette_versions <- c(
  V1_Original       = "#4F6B9A",   # blue
  V2_DecontX        = "#6B9BD2",   # light blue
  V3_DecontX_Micro  = "#2C4A7C"    # dark blue
)

# --- Palette: bulk datasets ---------------------------------------------------
palette_bulk <- c(
  Bulk_AB42 = "#35705A",   # green
  Bulk_BKO  = "#9A6A2F"    # amber
)

# --- Palette: direction (NES / log2FC) --------------------------------------
palette_direction <- c(
  Positive = "#A14E50",   # red
  Negative = "#4F6B9A"    # blue
)

# --- Palette: cell types ------------------------------------------------------
palette_celltypes <- c(
  Astrocyte       = "#35705A",
  Microglia       = "#9A6A2F",
  Neuron          = "#4F6B9A",
  Oligodendrocyte = "#7B5EA7",
  OPC             = "#A85D7E",
  Endothelial     = "#5A8A9A",
  Pericyte        = "#8A7A5A",
  Unassigned      = "#999999"
)

# --- Diverging palette for heatmaps ------------------------------------------
heatmap_palette <- function(n = 100) {
  colorRampPalette(c("navy", "white", "firebrick3"))(n)
}

# --- Set as default -----------------------------------------------------------
theme_set(theme_publication())
