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

| Script | What it does | Detail level |
|---|---|---|
| `00_shared_comparison_utils.R` | Reusable functions for rules 1–5: GO-ID join, NES-direction table, `core_enrichment` split/intersect, per-dataset log2FC bind, ±1.5-capped heatmap matrix. | **Full** — single source of truth for the rules. |
| `01_unified_GO_tally.R` | Tally significant GO-BP terms and functional groups for all four datasets under the same threshold; the overview counts (353/156/87; 911; 308). | Medium. |
| `02_AB42_vs_snRNA.R` | Direct-Aβ comparison: antigen-processing/MHC activation and reduced sterol/steroid biosynthesis — cross NES table and shared leading-edge genes for these modules. | **Full.** |
| `03_BKO_vs_AB42_vs_V3.R` | Three-way separation: Aβ42 ∩ BKO overlap (101 terms); BKO ∩ V3 but **not** Aβ42 (18 terms — 12 same-direction microtubule/cilium programs, 83 shared leading-edge genes / 62 reproduced; 6 opposite-direction). | **Full.** |
| `04_leading_edge_LFC_heatmap.R` | Generic leading-edge-gene × dataset log2FC heatmap built through `00` (direction only, capped ±1.5). | Medium/template. |

## Inputs / outputs

- **Input:** the DE + GSEA result tables produced by `01_snRNAseq/` and
  `02_bulkRNAseq/`.
- **Output:** cross-dataset GO/NES tables, shared leading-edge tables, and the matrices
  behind the comparison figures. Small final summaries go to `../results/tables/`.
