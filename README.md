# Manuscript figures: MF lines, larval tissues and wild-type development

Code that makes the data figures and supplementary tables for the manuscript "Embryos programmed to uniform cell fates define a small regulatory proteome operating long before specialization of the bulk proteome".
There is one Jupyter notebook per main figure. Code shared between notebooks is in `common.py`.

## Folder contents

```
manuscript_figures/
├── common.py                       shared paths, TMT normalisation, MF-line statistics, plotting helpers, tissue colours
├── fig2_MF_lines_3hpf.ipynb        Figure 2, Fig S2, Table S1
├── fig3_MF_timecourse.ipynb        Figure 3, Tables S2–S4
├── fig5_larval_tissues.ipynb       Figure 5, Fig S3, Tables S5–S6
├── fig6_CanS_development.ipynb     Figure 6, Figs S4–S5, Tables S7–S8
├── raw data/                       TMT protein quantification tables (inputs)
├── helper files/                   annotation and reference files (inputs)
└── outputs/                        everything the notebooks write (created on first run)
    ├── fig2/ fig3/ fig5/ fig6/     figure panels, each as .pdf (editable text) and .png (300 dpi)
    ├── supplementary_tables/       Tables S1–S8 as .csv
    └── intermediate/               files passed from fig3 and fig5 to fig6
```

### Input files

| File | Contents | Used by |
|---|---|---|
| `raw data/protein_quant_2228_cns15_b0.tsv` | Four MF lines, blastoderm (3 hpf), 4 replicates each | fig2, Table S1 |
| `raw data/protein_quant_2229_cns16_b4.tsv` | Four MF lines, blastoderm + 4 h (7 hpf), 4 replicates each | fig2 (Fig S2B), fig3 (Fig 3E, Table S4) |
| `raw data/protein_quant_2245_cns19b.tsv` | Four MF lines, time course at 1, 3, 5 and 7 hpf | fig3 |
| `raw data/protein_quant_2241_cns20_larva_tissues.tsv` | Larval brain, intestine, Malpighian tubules and imaginal discs | fig5 |
| `raw data/protein_quant_2242_cns21_devseries.tsv` | Wild-type (Canton-S) developmental series | fig6 |
| `helper files/genetable.xlsx` | Maps protein IDs (`ProteinId`) to gene symbols and FlyBase IDs | all |
| `helper files/gene_group_data_fb_2025_05.tsv.gz` | FlyBase gene groups (release FB2025_05) | fig2, fig5, fig6 |
| `helper files/figure3_go_reference.csv` | Saved GO terms for Fig 3C, used when g:Profiler can't be reached | fig3 |

## Running the notebooks

Run the notebooks in this order. **fig6 needs output from fig3 and fig5.**

```
fig2  →  fig3  →  fig5  →  fig6
```

| Notebook writes | Read by |
|---|---|
| `outputs/intermediate/mf_timecourse_clusters.csv` (fig3) | fig6 (Fig 6E) |
| `outputs/intermediate/larva_tissue_enriched_proteins.csv` (fig5) | fig6 (Figs 6B–E, S5) |

Open each notebook in Jupyter and choose *Run All*. You can also run them from the command line, inside this folder:

```bash
jupyter nbconvert --to notebook --execute --inplace fig2_MF_lines_3hpf.ipynb
jupyter nbconvert --to notebook --execute --inplace fig3_MF_timecourse.ipynb
jupyter nbconvert --to notebook --execute --inplace fig5_larval_tissues.ipynb
jupyter nbconvert --to notebook --execute --inplace fig6_CanS_development.ipynb
```

Start Jupyter from this folder (`manuscript_figures/`) so the notebooks can `import common`. All paths are set
relative to `common.py`, so the folder can be placed anywhere.

**Internet access:** Fig 3C and Fig 5E query the [g:Profiler](https://biit.cs.ut.ee/gprofiler/) web service for GO
enrichment. If g:Profiler can't be reached, Fig 3C falls back automatically to the saved terms in
`helper files/figure3_go_reference.csv`; Fig 5E has no offline fallback.

**Editing `common.py`:** if you change `common.py` while a notebook is open, restart the kernel to load the change.

### Software

The notebooks were run with Python 3.11 and these package versions:

| Package | Version |
|---|---|
| pandas | 2.3.3 |
| numpy | 2.3.3 |
| scipy | 1.16.2 |
| statsmodels | 0.14.6 |
| scikit-learn | 1.8.0 |
| matplotlib | 3.10.6 |
| seaborn | 0.13.2 |
| openpyxl | 3.1.5 (reads `genetable.xlsx`) |
| gprofiler-official | 1.0.0 |
| adjustText | 1.3.0 |

```bash
pip install pandas numpy scipy statsmodels scikit-learn matplotlib seaborn openpyxl gprofiler-official adjustText jupyter
```

## What each notebook makes

Schematic panels (2A, 3A left, 5A schematic and 6A) and Fig S1 (immunostaining) are not made by code.

### `fig2_MF_lines_3hpf.ipynb`: Figure 2

| Output | Description |
|---|---|
| Fig 2B | % of the proteome differentially abundant between each pair of lines, 3 hpf |
| Fig 2C | Kr, byn and zen across the four lines |
| Fig 2D | FlyBase gene-group violins |
| Fig S2A, S2B | Replicate concordance (best replicate pair per line) at 3 hpf and 7 hpf |
| Table S1 | All proteins at 3 hpf with log2FC, q-value and significance for the six comparisons |

### `fig3_MF_timecourse.ipynb`: Figure 3

| Output | Description |
|---|---|
| Fig 3A (right) | Ind, Byn and Peb profiles over time |
| Fig 3B | Shared K-means clusters (K = 12), median profile per line |
| Fig 3C | GO enrichment of clusters 3, 6, 8 and 9 |
| Fig 3D | Abundance of cluster 9 proteins at 7 hpf, per line |
| Fig 3E | % of the proteome differentially abundant between lines, 7 hpf |
| Table S2 | Time course, all proteins, with cluster |
| Table S3 | Proteins in the fate-specific clusters 6, 8 and 9 |
| Table S4 | All proteins at 7 hpf with log2FC, q-value and significance for the six comparisons |

### `fig5_larval_tissues.ipynb`: Figure 5

| Output | Description |
|---|---|
| Fig 5A (heatmap) | Median abundance of all proteins per tissue |
| Fig 5B | % of the proteome differing ≥ 2-fold between each pair of tissues |
| Fig 5C | FlyBase gene-group violins across tissues |
| Fig 5D | Tissue-prevalent (2–10×) and tissue-enriched (≥ 10×) proteins per tissue |
| Fig 5E | GO enrichment of tissue-enriched proteins |
| Fig S3 | Replicate concordance per tissue |
| Table S5 | All proteins with tissue medians, pairwise fold changes and adjusted p-values, and tissue calls |
| Table S6 | The 798 tissue-enriched proteins (subset of Table S5) |

### `fig6_CanS_development.ipynb`: Figure 6

| Output | Description |
|---|---|
| Fig 6B | K-means (K = 6) heatmap of the developmental time course |
| Fig 6C | Developmental onset of tissue-enriched proteins |
| Fig 6D | Brain-enriched proteins in the shared clusters |
| Fig 6E | Tissue-enriched proteins in MF time-course clusters 3, 6, 8 and 9 |
| Fig S4 | Ten FlyBase gene groups across development |
| Fig S5 | Shared clusters as in Fig 6D, one panel per tissue |
| Table S7 | Developmental time course with K-means cluster and tissue of enrichment |
| Table S8 | Tissue-enriched proteins detected in the time course, with shared cluster |

## Methods summary

- **TMT normalisation** (`common.normalize_tmt`). Steps, in order:
  1. Remove decoy and contaminant proteins.
  2. Keep proteins whose summed signal-to-noise across the channels used passes the threshold set in each notebook.
  3. Normalise each protein to its mean.
  4. Normalise each channel to its median.
  5. Normalise each protein to its mean again.
  6. Annotate with `genetable.xlsx`.
- **MF lines (Figs 2 and 3E, Tables S1 and S4).**
  - Welch t-test per protein, with Benjamini–Hochberg FDR.
  - log2 fold change from the replicate medians.
  - Significant: q < 0.05 and |log2FC| ≥ 1.
- **MF time course (Fig 3).**
  - K-means with K = 12 on the concatenated profiles of the four lines.
  - Clusters are renumbered 1–12 by size, largest first.
- **Larval tissues (Fig 5, Tables S5 and S6).**
  - Pairwise Welch t-tests between tissues, with BH-FDR across all pairs.
  - Zeros are replaced by half the smallest non-zero value.
  - **Tissue-enriched:** significantly higher than every other tissue, with a smallest log2FC ≥ 3.33 (≥ 10×).
  - **Tissue-prevalent:** significantly higher than every other tissue, with a smallest fold change of 2–10×.
  - In figures and tables the gut is called *intestine*. The intestine has three replicates (gutA excluded).
- **Wild-type development (Fig 6).**
  - K-means with K = 6 on standardised profiles of the whole time course (Fig 6B).
  - A separate K-means with K = 6 fit on all tissue-enriched proteins, shared between Fig 6D and Fig S5.
