# TREM2 PET scRNA-seq analysis

Start from processed Seurat objects. Existing clustering and UMAP coordinates are retained.

Set the working directory to this folder and edit paths in `00_settings.R`. Run `01_annotation.R` first, then run the other scripts as needed in R/RStudio.

- `00_settings.R`: shared settings and functions; loaded by each script.
- `01_annotation.R`: final cell annotations and annotated objects.
- `02_figures_and_proportions.R`: UMAPs, marker plots and cell proportions.
- `03_late_vs_early_DE_GO.R`: animal-level pseudobulk DEG, volcano plots and GO.
- `04_trem2_detection_DE_GO.R`: paired Trem2-positive/negative DEG and GO in DAM and Mono/Mac at each tumor stage.
- `05_gene_detection_percentages.R`: Trem2/Tfrc-positive percentages by sample, stage and cell type.

Annotation uses the original `clean.seurat`, `myeloid.seurat.clean` and `microglia.seurat.clean` objects in the R session, or their RDS files specified in `00_settings.R`. **Microglial subtype labels must precede renaming.** Final UMAP coordinates are reused from `all.seurat.clean` and `real.microglia.seurat.clean` when present, or optional coordinate files; otherwise input coordinates are retained.

Outputs are saved under `results/`. Rerunning annotation replaces its output objects. DEG/GO runs use separate output folders. Trem2 positivity means raw RNA counts >0. Check the summary CSVs for skipped or failed comparisons.
