#!/usr/bin/env Rscript

# SB28 glioblastoma single-cell RNA-seq analysis
#
# This script combines SAMPLE_DETECT.R and Untitled.R into one linear workflow.
# CD11b+ myeloid cells and GFP+ tumor cells were FACS-sorted, sequenced
# together, and analyzed jointly. The two populations separated clearly during
# initial clustering and were subsequently analyzed as separate subsets.


# =============================================================================
# 1. Packages and paths
# =============================================================================

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

# Change these paths before running the script.
cellranger_dir <- "PATH/TO/CellRanger_data"
output_dir <- "PATH/TO/output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "QC"), showWarnings = FALSE)
dir.create(file.path(output_dir, "batch"), showWarnings = FALSE)
dir.create(file.path(output_dir, "all_cells"), showWarnings = FALSE)
dir.create(file.path(output_dir, "microglia"), showWarnings = FALSE)
dir.create(file.path(output_dir, "microglia", "features"), showWarnings = FALSE)
dir.create(file.path(output_dir, "microglia", "DE"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tumor"), showWarnings = FALSE)


# =============================================================================
# 2. Read and merge Cell Ranger matrices
# =============================================================================

sample_names <- list.files(cellranger_dir)

all_mtx <- Read10X_h5(
  file.path(
    cellranger_dir,
    sample_names[1],
    "outs",
    "filtered_feature_bc_matrix.h5"
  )
)

colnames(all_mtx) <- paste(
  sample_names[1],
  colnames(all_mtx),
  sep = "_"
)

for (sample_name in sample_names[2:length(sample_names)]) {
  message("Loading ", sample_name)

  sample_mtx <- Read10X_h5(
    file.path(
      cellranger_dir,
      sample_name,
      "outs",
      "filtered_feature_bc_matrix.h5"
    )
  )

  colnames(sample_mtx) <- paste(
    sample_name,
    colnames(sample_mtx),
    sep = "_"
  )

  all_mtx <- cbind(all_mtx, sample_mtx)
}

all_seurat <- CreateSeuratObject(
  counts = all_mtx,
  project = "hospital"
)


# =============================================================================
# 3. Quality control
# =============================================================================

all_seurat[["percent.mt"]] <- PercentageFeatureSet(
  all_seurat,
  pattern = "^mt-"
)

all_seurat[["percent.ribo"]] <- PercentageFeatureSet(
  all_seurat,
  pattern = "^Rps|^Rpl"
)

pdf(
  file.path(output_dir, "QC", "QC_before_cutoff.pdf"),
  width = 5,
  height = 5
)

VlnPlot(
  all_seurat,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  pt.size = 0
)

dev.off()

all_seurat <- subset(
  all_seurat,
  subset = percent.mt < 10 &
    nCount_RNA <= 40000 &
    nFeature_RNA <= 6000
)

pdf(
  file.path(output_dir, "QC", "QC_after_cutoff.pdf"),
  width = 5,
  height = 5
)

VlnPlot(
  all_seurat,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  pt.size = 0
)

dev.off()

pdf(
  file.path(output_dir, "QC", "QC_per_sample.pdf"),
  width = 6,
  height = 4
)

VlnPlot(
  all_seurat,
  features = "nCount_RNA",
  group.by = "orig.ident",
  pt.size = 0
)

dev.off()


# =============================================================================
# 4. Initial joint analysis of all cells
# =============================================================================

# No integration or batch-correction method was used.
# The merged count matrix was normalized using SCTransform.

all_seurat <- SCTransform(
  all_seurat,
  verbose = TRUE,
  variable.features.rv.th = 1.4,
  return.only.var.genes = FALSE,
  variable.features.n = NULL
)

VariableFeatures(all_seurat) <- setdiff(
  VariableFeatures(all_seurat),
  str_subset(
    VariableFeatures(all_seurat),
    pattern = "^mt-"
  )
)

DefaultAssay(all_seurat) <- "SCT"

all_seurat <- RunPCA(
  all_seurat,
  npcs = 50,
  verbose = TRUE
)

pdf(
  file.path(output_dir, "all_cells", "elbow_all_cells.pdf"),
  width = 5,
  height = 5
)

ElbowPlot(all_seurat, ndims = 50)
dev.off()

all_seurat <- FindNeighbors(
  all_seurat,
  dims = 1:15
)

all_seurat <- FindClusters(
  all_seurat,
  resolution = 0.3
)

all_seurat <- RunUMAP(
  all_seurat,
  dims = 1:15,
  do.fast = TRUE,
  n.neighbors = 50
)

p <- DimPlot(
  all_seurat,
  label = TRUE
) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )

ggsave(
  file.path(output_dir, "all_cells", "UMAP_all_cells.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)


# =============================================================================
# 5. Per-sample UMAP overlays
# =============================================================================

umap_metadata <- cbind(
  as.data.frame(all_seurat@meta.data),
  as.data.frame(Embeddings(all_seurat, reduction = "umap"))
)

for (sample_name in unique(umap_metadata$orig.ident)) {
  p <- ggplot() +
    geom_point(
      data = umap_metadata,
      aes(x = UMAP_1, y = UMAP_2),
      color = "gray88",
      size = 0.5
    ) +
    geom_point(
      data = umap_metadata[
        umap_metadata$orig.ident == sample_name,
      ],
      aes(x = UMAP_1, y = UMAP_2),
      color = "gray6",
      size = 0.5
    ) +
    ggtitle(sample_name) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )

  ggsave(
    file.path(
      output_dir,
      "batch",
      paste0(sample_name, "_UMAP.png")
    ),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

all_cell_markers <- FindAllMarkers(all_seurat)

write.table(
  all_cell_markers,
  file.path(output_dir, "all_cells", "markers_all_cells.txt"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

saveRDS(
  all_seurat,
  file.path(output_dir, "all_cells", "all_seurat_initial.rds")
)


# =============================================================================
# 6. Select samples used in the clean analysis
# =============================================================================

# The original clean analysis used two samples per condition.

clean_seurat <- subset(
  all_seurat,
  subset = orig.ident %in% c(
    "E1",
    "E2",
    "L1",
    "L2",
    "S1",
    "S2"
  )
)

clean_counts <- GetAssayData(
  clean_seurat,
  assay = "RNA",
  slot = "counts"
)

clean_metadata <- clean_seurat@meta.data

clean_seurat <- CreateSeuratObject(
  counts = clean_counts,
  meta.data = clean_metadata
)

clean_seurat <- SCTransform(
  clean_seurat,
  verbose = TRUE,
  variable.features.rv.th = 1.4,
  return.only.var.genes = FALSE,
  variable.features.n = NULL
)

VariableFeatures(clean_seurat) <- setdiff(
  VariableFeatures(clean_seurat),
  str_subset(
    VariableFeatures(clean_seurat),
    pattern = "^mt-"
  )
)

clean_seurat <- RunPCA(
  clean_seurat,
  npcs = 50,
  verbose = TRUE
)

pdf(
  file.path(output_dir, "all_cells", "elbow_clean_cells.pdf"),
  width = 5,
  height = 5
)

ElbowPlot(clean_seurat, ndims = 50)
dev.off()

clean_seurat <- FindNeighbors(
  clean_seurat,
  dims = 1:23
)

clean_seurat <- FindClusters(
  clean_seurat,
  resolution = 0.1
)

clean_seurat <- RunUMAP(
  clean_seurat,
  dims = 1:23,
  do.fast = TRUE,
  n.neighbors = 50
)


# =============================================================================
# 7. Annotate major cell populations
# =============================================================================

# These cluster mappings are specific to the original dataset.

clean_seurat$cluster_number <- as.numeric(
  as.character(clean_seurat$seurat_clusters)
)

clean_seurat$major <- "Mac/microglia"
clean_seurat$major[
  clean_seurat$cluster_number == 0
] <- "Tumor"
clean_seurat$major[
  clean_seurat$cluster_number == 5
] <- "NK"
clean_seurat$major[
  clean_seurat$cluster_number == 6
] <- "Monocyte"

clean_seurat$condition <- "Control"
clean_seurat$condition[
  clean_seurat$orig.ident %in% c("E1", "E2")
] <- "Early_stage"
clean_seurat$condition[
  clean_seurat$orig.ident %in% c("L1", "L2")
] <- "Late_stage"

clean_seurat$condition <- factor(
  clean_seurat$condition,
  levels = c(
    "Control",
    "Early_stage",
    "Late_stage"
  )
)

Idents(clean_seurat) <- "major"

major_colors <- c(
  "Mac/microglia" = "#2F7FC1",
  "Tumor" = "#D8383A",
  "NK" = "#96C37D",
  "Monocyte" = "#F3D266"
)

p <- DimPlot(
  clean_seurat,
  group.by = "major",
  cols = major_colors
) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

ggsave(
  file.path(output_dir, "all_cells", "UMAP_major_cell_types.pdf"),
  p,
  width = 8,
  height = 5
)

p <- ggplot(
  clean_seurat@meta.data,
  aes(x = condition, fill = major)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(
    labels = percent_format(),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = major_colors) +
  xlab("") +
  ylab("Fraction of cells") +
  theme_classic()

ggsave(
  file.path(output_dir, "all_cells", "major_cell_fraction.pdf"),
  p,
  width = 4.5,
  height = 5
)

saveRDS(
  clean_seurat,
  file.path(output_dir, "all_cells", "clean_seurat.rds")
)


# =============================================================================
# 8. Mac/microglia subclustering
# =============================================================================

microglia_seurat <- subset(
  clean_seurat,
  subset = major == "Mac/microglia"
)

microglia_counts <- GetAssayData(
  microglia_seurat,
  assay = "RNA",
  slot = "counts"
)

microglia_metadata <- microglia_seurat@meta.data

microglia_seurat <- CreateSeuratObject(
  counts = microglia_counts,
  meta.data = microglia_metadata
)

microglia_seurat <- SCTransform(
  microglia_seurat,
  verbose = TRUE,
  variable.features.rv.th = 1.4,
  return.only.var.genes = FALSE,
  variable.features.n = NULL
)

VariableFeatures(microglia_seurat) <- setdiff(
  VariableFeatures(microglia_seurat),
  str_subset(
    VariableFeatures(microglia_seurat),
    pattern = "^mt-"
  )
)

microglia_seurat <- RunPCA(
  microglia_seurat,
  npcs = 50,
  verbose = TRUE
)

microglia_seurat <- FindNeighbors(
  microglia_seurat,
  dims = 1:20
)

microglia_seurat <- FindClusters(
  microglia_seurat,
  resolution = 0.2
)

microglia_seurat <- RunUMAP(
  microglia_seurat,
  dims = 1:20,
  do.fast = TRUE,
  n.neighbors = 50
)

p <- DimPlot(
  microglia_seurat,
  label = TRUE
) +
  theme_classic()

ggsave(
  file.path(output_dir, "microglia", "UMAP_microglia_initial.png"),
  p,
  width = 6,
  height = 4,
  dpi = 300
)

microglia_markers_initial <- FindAllMarkers(microglia_seurat)

write.table(
  microglia_markers_initial,
  file.path(
    output_dir,
    "microglia",
    "markers_microglia_initial.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# =============================================================================
# 9. Final Mac/microglia clustering
# =============================================================================



microglia_counts <- GetAssayData(
  microglia_seurat,
  assay = "RNA",
  slot = "counts"
)

microglia_metadata <- microglia_seurat@meta.data

microglia_clean <- CreateSeuratObject(
  counts = microglia_counts,
  meta.data = microglia_metadata
)

microglia_clean <- SCTransform(
  microglia_clean,
  verbose = TRUE,
  variable.features.rv.th = 1.4,
  return.only.var.genes = FALSE,
  variable.features.n = NULL
)

VariableFeatures(microglia_clean) <- setdiff(
  VariableFeatures(microglia_clean),
  str_subset(
    VariableFeatures(microglia_clean),
    pattern = "^mt-"
  )
)

microglia_clean <- RunPCA(
  microglia_clean,
  npcs = 50,
  verbose = TRUE
)

microglia_clean <- FindNeighbors(
  microglia_clean,
  dims = 1:15
)

microglia_clean <- FindClusters(
  microglia_clean,
  resolution = 0.5
)

microglia_clean <- RunUMAP(
  microglia_clean,
  dims = 1:15,
  do.fast = TRUE
)


# =============================================================================
# 10. Annotate Homeostatic, pre-DAM, and DAM clusters
# =============================================================================

# These cluster mappings are specific to the original analysis.

microglia_clean$cluster_number <- as.numeric(
  as.character(microglia_clean$seurat_clusters)
)

microglia_clean$subtype <- "Homeostatic"

microglia_clean$subtype[
  microglia_clean$cluster_number %in% c(2, 12)
] <- "pre-DAM"

microglia_clean$subtype[
  microglia_clean$cluster_number %in% c(0, 8)
] <- "DAM1"

microglia_clean$subtype[
  microglia_clean$cluster_number == 5
] <- "DAM2"

microglia_clean$subtype[
  microglia_clean$cluster_number == 6
] <- "DAM3"

microglia_clean$subtype[
  microglia_clean$cluster_number == 7
] <- "DAM4"

microglia_clean$subtype[
  microglia_clean$cluster_number == 1
] <- "DAM5"

microglia_clean$subtype[
  microglia_clean$cluster_number == 11
] <- "DAM6"

microglia_clean$subtype <- factor(
  microglia_clean$subtype,
  levels = c(
    "Homeostatic",
    "pre-DAM",
    "DAM1",
    "DAM2",
    "DAM3",
    "DAM4",
    "DAM5",
    "DAM6"
  )
)

Idents(microglia_clean) <- "subtype"

microglia_colors <- c(
  "Homeostatic" = "#4E79A7",
  "pre-DAM" = "#59A14F",
  "DAM1" = "#F28E2B",
  "DAM2" = "#E15759",
  "DAM3" = "#76B7B2",
  "DAM4" = "#EDC948",
  "DAM5" = "#B07AA1",
  "DAM6" = "#FF9DA7"
)

p <- DimPlot(
  microglia_clean,
  group.by = "subtype",
  cols = microglia_colors,
  label = TRUE
) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  )

ggsave(
  file.path(output_dir, "microglia", "UMAP_DAM_subtypes.png"),
  p,
  width = 6,
  height = 5,
  dpi = 300
)

p <- DimPlot(
  microglia_clean,
  group.by = "subtype",
  split.by = "condition",
  cols = microglia_colors
) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

ggsave(
  file.path(
    output_dir,
    "microglia",
    "UMAP_DAM_split_by_condition.png"
  ),
  p,
  width = 11,
  height = 4,
  dpi = 300
)

p <- ggplot(
  microglia_clean@meta.data,
  aes(x = condition, fill = subtype)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(
    labels = percent_format(),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = microglia_colors) +
  xlab("") +
  ylab("Fraction of Mac/microglia") +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "microglia",
    "DAM_fraction_by_condition.pdf"
  ),
  p,
  width = 4.5,
  height = 5
)


# =============================================================================
# 11. DAM marker analysis
# =============================================================================

microglia_markers <- FindAllMarkers(
  microglia_clean,
  min.pct = 0.1
)

write.table(
  microglia_markers,
  file.path(
    output_dir,
    "microglia",
    "markers_DAM_subtypes.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

top_markers <- microglia_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  )

top_genes <- unique(top_markers$gene)

average_expression <- AverageExpression(
  microglia_clean,
  assays = "SCT",
  features = top_genes,
  return.seurat = TRUE
)

heatmap_matrix <- GetAssayData(
  average_expression,
  assay = "SCT",
  slot = "scale.data"
)

pdf(
  file.path(
    output_dir,
    "microglia",
    "heatmap_DAM_markers.pdf"
  ),
  width = 4,
  height = 8
)

pheatmap(
  heatmap_matrix,
  breaks = seq(-2, 2, 0.1),
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(40),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA
)

dev.off()


# =============================================================================
# 12. Key gene expression plots
# =============================================================================

key_genes <- c(
  "Trem2",
  "Tspo",
  "P2ry12",
  "Apoe",
  "Spp1",
  "Cd74",
  "Arg1",
  "H2-Aa",
  "H2-Ab1"
)

for (gene in key_genes) {
  p <- FeaturePlot(
    microglia_clean,
    features = gene,
    cols = brewer.pal(5, "OrRd"),
    raster = FALSE
  ) +
    ggtitle(gene) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )

  ggsave(
    file.path(
      output_dir,
      "microglia",
      "features",
      paste0(gene, "_FeaturePlot.pdf")
    ),
    p,
    width = 6,
    height = 5
  )
}

# Trem2 is interpreted as a myeloid-enriched and state-dependent marker,
# rather than a marker exclusive to DAM clusters.

p <- VlnPlot(
  microglia_clean,
  features = "Trem2",
  group.by = "subtype",
  pt.size = 0
) +
  xlab("") +
  ylab("Trem2 expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(
    output_dir,
    "microglia",
    "Trem2_expression_by_subtype.pdf"
  ),
  p,
  width = 6,
  height = 4
)

p <- DotPlot(
  microglia_clean,
  features = "Trem2",
  group.by = "subtype"
) +
  RotatedAxis() +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "microglia",
    "Trem2_dotplot_by_subtype.pdf"
  ),
  p,
  width = 5,
  height = 3
)

saveRDS(
  microglia_clean,
  file.path(
    output_dir,
    "microglia",
    "microglia_DAM_annotated.rds"
  )
)


# =============================================================================
# 13. Differential expression: Homeostatic versus DAM clusters
# =============================================================================

microglia_de_counts <- GetAssayData(
  microglia_clean,
  assay = "SCT",
  slot = "counts"
)

microglia_de <- CreateSeuratObject(
  counts = microglia_de_counts,
  meta.data = microglia_clean@meta.data
)

microglia_de <- NormalizeData(
  microglia_de,
  normalization.method = "LogNormalize"
)

Idents(microglia_de) <- "subtype"

comparison_clusters <- c(
  "pre-DAM",
  "DAM1",
  "DAM2",
  "DAM3",
  "DAM4",
  "DAM5",
  "DAM6"
)

for (comparison_cluster in comparison_clusters) {
  de_table <- FindMarkers(
    microglia_de,
    ident.1 = "Homeostatic",
    ident.2 = comparison_cluster,
    min.pct = 0,
    test.use = "DESeq2",
    logfc.threshold = 0,
    max.cells.per.ident = 50
  )

  de_table$gene <- rownames(de_table)

  write.table(
    de_table,
    file.path(
      output_dir,
      "microglia",
      "DE",
      paste0(
        "Homeostatic_vs_",
        comparison_cluster,
        ".txt"
      )
    ),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
}


# =============================================================================
# 14. Differential expression: early DAM1 versus late DAM clusters
# =============================================================================

microglia_de$subtype_condition <- paste0(
  microglia_de$subtype,
  "_",
  microglia_de$condition
)

Idents(microglia_de) <- "subtype_condition"

late_DAM_clusters <- c(
  "DAM1_Late_stage",
  "DAM2_Late_stage",
  "DAM3_Late_stage",
  "DAM4_Late_stage",
  "DAM5_Late_stage",
  "DAM6_Late_stage"
)

for (late_DAM_cluster in late_DAM_clusters) {
  de_table <- FindMarkers(
    microglia_de,
    ident.1 = "DAM1_Early_stage",
    ident.2 = late_DAM_cluster,
    min.pct = 0,
    test.use = "DESeq2",
    logfc.threshold = 0,
    max.cells.per.ident = 50
  )

  de_table$gene <- rownames(de_table)

  write.table(
    de_table,
    file.path(
      output_dir,
      "microglia",
      "DE",
      paste0(
        "DAM1_Early_vs_",
        late_DAM_cluster,
        ".txt"
      )
    ),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
}


# =============================================================================
# 15. Tumor-cell subclustering
# =============================================================================

tumor_seurat <- subset(
  clean_seurat,
  subset = major == "Tumor"
)

tumor_counts <- GetAssayData(
  tumor_seurat,
  assay = "RNA",
  slot = "counts"
)

tumor_metadata <- tumor_seurat@meta.data

tumor_seurat <- CreateSeuratObject(
  counts = tumor_counts,
  meta.data = tumor_metadata
)

tumor_seurat <- SCTransform(
  tumor_seurat,
  verbose = TRUE,
  variable.features.rv.th = 1.4,
  return.only.var.genes = FALSE,
  variable.features.n = NULL
)

VariableFeatures(tumor_seurat) <- setdiff(
  VariableFeatures(tumor_seurat),
  str_subset(
    VariableFeatures(tumor_seurat),
    pattern = "^mt-"
  )
)

tumor_seurat <- RunPCA(
  tumor_seurat,
  npcs = 50,
  verbose = TRUE
)

tumor_seurat <- FindNeighbors(
  tumor_seurat,
  dims = 1:15
)

tumor_seurat <- FindClusters(
  tumor_seurat,
  resolution = 0.1
)

tumor_seurat <- RunUMAP(
  tumor_seurat,
  dims = 1:15,
  do.fast = TRUE,
  n.neighbors = 50
)

tumor_seurat$subtype <- paste0(
  "tumor",
  as.numeric(
    as.character(tumor_seurat$seurat_clusters)
  ) + 1
)

Idents(tumor_seurat) <- "subtype"

p <- DimPlot(
  tumor_seurat,
  group.by = "subtype",
  label = TRUE
) +
  theme_classic()

ggsave(
  file.path(output_dir, "tumor", "UMAP_tumor_subtypes.png"),
  p,
  width = 6,
  height = 4,
  dpi = 300
)

p <- DimPlot(
  tumor_seurat,
  group.by = "subtype",
  split.by = "condition",
  label = TRUE
) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "tumor",
    "UMAP_tumor_split_by_condition.png"
  ),
  p,
  width = 10,
  height = 4,
  dpi = 300
)

tumor_markers <- FindAllMarkers(tumor_seurat)

write.table(
  tumor_markers,
  file.path(
    output_dir,
    "tumor",
    "markers_tumor_subtypes.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# =============================================================================
# 16. Tspo and Tfrc expression
# =============================================================================

p <- FeaturePlot(
  tumor_seurat,
  features = "Tspo",
  cols = brewer.pal(5, "OrRd"),
  raster = FALSE
) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "tumor",
    "Tspo_tumor_FeaturePlot.pdf"
  ),
  p,
  width = 6,
  height = 5
)

p <- FeaturePlot(
  clean_seurat,
  features = "Tfrc",
  cols = brewer.pal(5, "OrRd"),
  raster = FALSE
) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "all_cells",
    "Tfrc_FeaturePlot.pdf"
  ),
  p,
  width = 6,
  height = 5
)

p <- DotPlot(
  clean_seurat,
  features = "Tfrc",
  group.by = "major"
) +
  RotatedAxis() +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "all_cells",
    "Tfrc_dotplot_major_cells.pdf"
  ),
  p,
  width = 5,
  height = 3
)


# =============================================================================
# 17. Save final files
# =============================================================================

saveRDS(
  tumor_seurat,
  file.path(
    output_dir,
    "tumor",
    "tumor_subclustered.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

message("Analysis complete.")
