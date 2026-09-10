# 03_cross_dataset — cross-modality comparison

This folder formalises **how the snRNA versions and the two bulk experiments are put
on the same footing**. It is written to be read carefully: the matching rules below are
what make the comparison valid.

## Comparison rules (implemented in `00_shared_comparison_utils.R`)

1. **Terms are matched by GO ID** (not by term text), under one identical threshold
   (`padj<0.05 & |NES|≥1.5`) and one 10-group dictionary.
2. **Direction is compared on NES** (sign only): same-direction vs opposite-direction.
3. **Shared leading-edge genes** = intersection of the two sides' `core_enrichment`
   strings (split on `/`) for the same GO term.
4. **Gene-level log2FC stays modality-specific** — DESeq2 log2FC for bulk, FindMarkers
   avg_log2FC for snRNA. They are used for **direction only**, never compared as
   effect sizes across platforms.
5. Combined leading-edge heatmaps use a fixed colour cap of **±1.5** for visual parity.

## Inputs / outputs

- **Input:** the DE + GSEA result tables produced by `01_snRNAseq/` and
  `02_bulkRNAseq/`.
- **Output:** cross-dataset GO/NES tables, shared leading-edge tables, and the matrices
  behind the comparison figures. Small final summaries go to `../results/tables/`.
