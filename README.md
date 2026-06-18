# Single-Cell RNA-seq Analysis of the SB28 Glioblastoma Microenvironment

This repository contains the R analysis workflow used to characterize tumor-associated myeloid cells and tumor cells during SB28 glioblastoma progression.

The analysis accompanies the study:

> **Transport vehicle-mediated TREM2 PET maps myeloid cells in the glioblastoma microenvironment**

## Data availability

The single-cell RNA-seq dataset is available from the NCBI Gene Expression Omnibus:

**GEO accession: [GSE317439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317439)**

Replace `GSEXXXX` with the final GEO accession number before publishing this repository.

## Experimental design

CD11b-positive TAMs/microglia and GFP-positive SB28 tumor cells were isolated by FACS. The sorted populations were sequenced together and analyzed jointly during the initial Seurat workflow.

Tumor and myeloid populations separated clearly during initial clustering. Mac/microglia and tumor cells were subsequently subsetted and reclustered separately.

The final clean analysis included:

| Condition | Samples |
|---|---|
| Sham/control | S1, S2 |
| Early-stage tumor | E1, E2 |
| Late-stage tumor | L1, L2 |

## Analysis script

The main analysis is provided in:

[`scRNAseq_analysis_readable.R`](scRNAseq_analysis_readable.R)

The script follows a linear workflow:

1. Read and merge Cell Ranger filtered count matrices.
2. Calculate mitochondrial and ribosomal read percentages.
3. Apply quality-control filtering.
4. Normalize the merged dataset using SCTransform.
5. Perform PCA, graph-based clustering, and UMAP.
6. Annotate major cell populations.
7. Subset and recluster Mac/microglia.
8. Annotate Homeostatic, pre-DAM, and DAM1-DAM6 populations.
9. Generate marker, composition, and gene-expression plots.
10. Perform differential expression analyses.
11. Subset and recluster tumor cells.

## Requirements

The analysis was performed in R using the following packages:

```r
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggrepel)
library(scales)
```

Cell Ranger output for each sample must contain:

```text
SAMPLE_ID/
└── outs/
    └── filtered_feature_bc_matrix.h5
```

## Usage

Edit the input and output paths near the beginning of the script:

```r
cellranger_dir <- "PATH/TO/CellRanger_data"
output_dir <- "PATH/TO/output"
```

Run the analysis from the command line:

```bash
Rscript scRNAseq_analysis_readable.R
```

## Main analysis parameters

### Quality control

- Mitochondrial reads: `<10%`
- Total UMI count: `<=40,000`
- Detected genes: `<=6,000`

### Initial joint analysis

- Normalization: SCTransform
- Variable-feature residual variance threshold: `1.4`
- PCA components calculated: `50`
- Dimensions used for neighbors and UMAP: `1:15`
- Clustering resolution: `0.3`
- UMAP neighbors: `50`

### Clean all-cell analysis

- Dimensions used for neighbors and UMAP: `1:23`
- Clustering resolution: `0.1`
- UMAP neighbors: `50`

### Mac/microglia analysis

- Initial subclustering dimensions: `1:20`
- Initial clustering resolution: `0.2`
- Final clustering dimensions: `1:15`
- Final clustering resolution: `0.5`

No Seurat integration or external batch-correction method was applied. Samples were merged before SCTransform normalization, and per-sample UMAP overlays were generated to inspect sample distributions.

## Outputs

The script creates separate output directories for:

- quality-control plots;
- all-cell UMAP and marker results;
- per-sample UMAP overlays;
- Mac/microglia and DAM analyses;
- differential expression results;
- tumor-cell analyses;
- serialized Seurat objects;
- R session information.

## Important annotation note

Cluster-to-cell-type mappings are specific to this dataset and should be checked against canonical marker expression when the analysis is rerun.

The original exploratory analysis also included interactive `CellSelector`-based removal of selected cells. Because an interactive selection cannot be reconstructed from code alone, the selected cell barcodes or the corresponding final Seurat object should be provided if exact reproduction of those manually curated figures is required.

## Interpretation of TREM2

Trem2 expression is observed across multiple myeloid states, including Homeostatic/pre-DAM populations and DAM clusters. TREM2 should therefore be interpreted as a myeloid-enriched and state-dependent marker rather than a marker exclusive to one DAM population.

## Citation

Please cite the associated manuscript and GEO dataset when using this code or data. Full citation information will be added following publication.
