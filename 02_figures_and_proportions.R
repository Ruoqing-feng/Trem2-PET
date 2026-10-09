source("00_settings.R")
out = new_dir("figures")
all.cells = read_object("all.cells")
myeloid = read_object("myeloid")
microglia = read_object("microglia")
clean_umap = theme(axis.title = element_blank(), axis.text = element_blank(),
                   axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank())
save_plot = function(p, name, w = 7, h = 6) ggsave(file.path(out, name), p, width = w, height = h)
save_plot(DimPlot(all.cells, group.by = "major", cols = c(Myeloid = "#2F7FC1", Tumor = "#D8383A", NK = "#96C37D")) + clean_umap,
          "All_cells_major_UMAP.pdf")
save_plot(DimPlot(all.cells, group.by = "myeloid_subtype2", cols = celltype.colors) + clean_umap,
          "All_cells_annotations_UMAP.pdf", 10, 6)
save_plot(DimPlot(all.cells, group.by = "myeloid_subtype2", split.by = "stage", cols = celltype.colors, ncol = 3) + clean_umap,
          "All_cells_split_stage_UMAP.pdf", 18, 6)
save_plot(DimPlot(all.cells, group.by = "myeloid_subtype2", split.by = "orig.ident", cols = celltype.colors, ncol = 3) + clean_umap,
          "All_cells_split_sample_UMAP.pdf", 18, 10)
save_plot(DimPlot(microglia, group.by = "subtype", cols = micro.colors, label = TRUE) + clean_umap,
          "Microglial_subclusters_UMAP.pdf")
save_plot(DimPlot(microglia, group.by = "subtype", split.by = "stage", cols = micro.colors, ncol = 3) + clean_umap,
          "Microglial_subclusters_split_stage_UMAP.pdf", 15, 5)

# Composition denominator: recovered myeloid cells in each sample (not all sorted cells).
composition = function(x, label) {
  md = x[[]]
  types = sort(unique(as.character(md[[label]])))
  counts = as.data.frame(table(sample = factor(md$orig.ident, levels = samples$sample),
                               celltype = factor(md[[label]], levels = types)))
  names(counts)[3] = "n_cells"
  counts$stage = factor(unname(sample.stage[as.character(counts$sample)]), levels = c("Sham", "Early", "Late"))
  per.sample = counts %>% group_by(stage, sample) %>%
    mutate(total_cells = sum(n_cells), percentage = ifelse(total_cells > 0, 100 * n_cells / total_cells, NA_real_)) %>% ungroup()
  per.stage = per.sample %>% group_by(stage, celltype) %>%
    summarise(n_cells = sum(n_cells), total_cells = sum(total_cells),
      percent_pooled = ifelse(total_cells > 0, 100 * n_cells / total_cells, NA_real_),
      percent_sample_mean = if (all(is.na(percentage))) NA_real_ else mean(percentage, na.rm = TRUE), .groups = "drop")
  list(sample = per.sample, stage = per.stage)
}
for (kind in c("myeloid", "microglia")) {
  x = get(kind)
  comp = composition(x, if (kind == "myeloid") "sub4" else "subtype")
  save_table(comp$sample, file.path(out, paste0(kind, "_composition_by_sample.csv")))
  save_table(comp$stage, file.path(out, paste0(kind, "_composition_by_stage.csv")))
  p = ggplot(comp$sample, aes(sample, percentage / 100, fill = celltype)) +
    geom_col() + facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    scale_fill_manual(values = if (kind == "myeloid") celltype.colors else micro.colors) +
    theme_classic() + labs(x = NULL, y = "Cell proportion", fill = NULL)
  save_plot(p, paste0(kind, "_composition.pdf"), 10, 5)
}

marker_heatmap = function(x, label, genes, order, limit, file) {
  if (!"SCT" %in% SeuratObject::Assays(x)) stop("SCT assay required for marker heatmaps.")
  DefaultAssay(x) = "SCT"
  Idents(x) = label
  genes.present = intersect(genes, rownames(x[["SCT"]]))
  if (length(setdiff(genes, genes.present))) warning("Markers missing from SCT: ", paste(setdiff(genes, genes.present), collapse = ", "))
  if (!length(genes.present)) stop("No requested heatmap markers in SCT.")
  avg = AverageExpression(x, assays = "SCT", features = genes.present, return.seurat = TRUE)
  scaled = SeuratObject::LayerData(avg[["SCT"]], layer = "scale.data")
  # Seurat may sanitize group names in AverageExpression column names.
  wanted = gsub("_", "-", order)
  cols = match(wanted, colnames(scaled))
  if (anyNA(cols)) {
    cols = match(order, colnames(scaled))
  }
  cols = cols[!is.na(cols)]
  if (!length(cols)) stop("Could not match heatmap groups.")
  scaled = scaled[genes.present, cols, drop = FALSE]
  write.csv(scaled, file.path(out, paste0(file, "_values.csv")))
  breaks = seq(-limit, limit, length.out = 61)
  p = pheatmap::pheatmap(scaled, breaks = breaks,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
    cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, silent = TRUE)
  save_plot(p$gtable, paste0(file, ".pdf"), 5, max(5, 0.22 * nrow(scaled) + 2))
  if (cfg$make.markers) {
    markers = FindAllMarkers(x, assay = "SCT", min.pct = 0.1, only.pos = TRUE, test.use = "wilcox")
    save_table(markers, file.path(out, paste0(file, "_FindAllMarkers.csv")))
  }
}
marker_heatmap(myeloid, "sub4", myeloid.genes, celltype.order, 2, "Myeloid_markers")
marker_heatmap(microglia, "subtype", micro.genes, micro.order, 1.5, "Microglial_markers")

DefaultAssay(microglia) = "SCT"
for (gene in intersect(unique(c(micro.genes, expression.genes)), rownames(microglia[["SCT"]]))) {
  save_plot(FeaturePlot(microglia, features = gene, reduction = "umap", cols = c("grey90", "firebrick3")) + clean_umap,
            paste0("Microglial_feature_", safe_name(gene), ".pdf"))
}
DefaultAssay(all.cells) = "SCT"
for (gene in intersect(unique(c(myeloid.genes, micro.genes, expression.genes)), rownames(all.cells[["SCT"]]))) {
  save_plot(FeaturePlot(all.cells, features = gene, reduction = "umap", cols = c("grey90", "firebrick3")) + clean_umap,
            paste0("Feature_", safe_name(gene), ".pdf"))
}
for (gene in expression.genes) {
  if (!gene %in% rownames(all.cells[["SCT"]])) next
  save_plot(VlnPlot(all.cells, features = gene, group.by = "myeloid_subtype2", pt.size = 0,
                   cols = celltype.colors) + RotatedAxis(), paste0("Violin_", gene, ".pdf"), 11, 5)
}
all.cells$major_stage = paste(all.cells$major, all.cells$stage, sep = "_")
Idents(all.cells) = "major_stage"
present = intersect(expression.genes, rownames(all.cells[["SCT"]]))
if (length(present)) save_plot(DotPlot(all.cells, features = present), "Trem2_Tfrc_dotplot.pdf", 6, 5)
writeLines("Feature/violin/dot plots use SCT expression; raw-RNA detection percentages are produced by script 05. DotPlot colors are scaled per gene and do not establish cross-gene fold ratios.", file.path(out, "expression_plot_scale.txt"))
finish_run(out)
