# 04_figures — reusable plotting templates only

This folder deliberately contains **templates, not the dozens of one-off historical
plot scripts**. Each file exposes a reusable function (or a clearly marked template
block) with the expected input structure documented, so any result table from
`01`–`03` can be plotted consistently. Final figure assembly for the manuscript stays
outside the repository; these templates guarantee a single theme/palette.

| File | Template for |
|---|---|
| `theme_palette.R` | Unified `ggplot2` theme and the fixed palette: three version colours (V1/V2/V3), two bulk colours (Aβ42 green, BKO amber), direction colours, font sizes. **Loaded by every other template.** |
| `template_umap_dotplot.R` | Cell-type UMAP and marker dot/feature plots (snRNA annotation QC). |
| `template_GSEA_nes_bar.R` | Horizontal NES bar charts and the 10-functional-group stacked bars. |
| `template_leading_edge_heatmap.R` | `pheatmap` leading-edge × dataset log2FC heatmap (fixed ±1.5 cap, consistent annotation). |
| `template_overlap_venn_upSet.R` | GO-ID overlap / UpSet-style summary across datasets and versions. |

## Convention

- Every template runs on a **small documented example / expected column names**, so it
  can be sourced without the full data; replace the example with a result table.
- No paths to personal drives; output goes to `../results/figures/` (git-ignored).
