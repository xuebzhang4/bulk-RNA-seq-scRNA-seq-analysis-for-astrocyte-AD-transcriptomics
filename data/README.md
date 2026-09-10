# Data

Only small **metadata/sample sheets** are tracked in this repository. Raw counts and
large intermediate objects (`.rds`, `.h5`, count matrices) are **not** committed —
download them into `data/raw/` (git-ignored) following the table below, or retrieve the
frozen archive from Zenodo (DOI to be added).

## Datasets used

| Label | Modality | Contrast | Samples | Source |
|---|---|---|---|---|
| snRNA | single-nucleus RNA-seq, mouse cortex | WT vs 5xFAD (saline arms only) | WT-saline `GSM7092584`; 5xFAD-saline `GSM7092586` | **GEO GSE227157** (Fatmi *et al.*, *Aging* 2024;16:3137–3159, PMID 38385967). The two APC arms (`GSM7092585`, `GSM7092587`) are excluded. |
| bulk Aβ42 | bulk RNA-seq, purified primary astrocytes | Aβ42 48 h vs untreated | `A_Beta_1–3` vs `Astro_1–3` (n=3 vs 3) | In-house. **GEO accession: to be added (not yet deposited).** |
| bulk BKO | bulk RNA-seq, purified astrocytes | Bmal1-KO vs WT | `Astro_BKO_1–3` (`GSM9610294–296`) vs `Astro_WT_1–4` (`GSM9610297–300`), n=3 vs 4 | **GEO GSE325658** (in-house prior project). |

Sample sheets with the exact group labels used by the scripts are in `metadata/`:

- `metadata/snRNA_samples.csv`
- `metadata/bulk_AB42_samples.csv`
- `metadata/bulk_BKO_samples.csv`

## Expected local layout (git-ignored)

```
data/raw/
  gse227157/                 # 10x outputs per sample (barcodes/features/matrix)
  bulk_AB42/Astrocyte_Only_Clean_Counts.csv
  bulk_BKO/                  # count matrix from GSE325658 (Astro_BKO / Astro_WT)
```

> Note: GSE206081 is a separate **neuron** project (early/late Aβ, excitotoxicity,
> neuronal BKO/WT) and is **not** used here.
