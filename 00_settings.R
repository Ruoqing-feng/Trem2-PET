# Shared settings and functions. Each analysis script sources this file.
# Set the working directory to this folder; edit input/output paths here.
cfg = list(
  input.dir = "data/inputs",
  output.dir = "results",
  inputs = c(all = "clean_seurat.rds", myeloid = "myeloid_clean.rds",
             microglia = "microglia_clean.rds"),
  all.umap = "all_umap.rds", microglia.umap = "microglia_umap.rds",
  stage.clusters = c("DAM", "Mono/Mac", "Homeostatic/early activated microglia", "MHC-II+ DC"),
  trem2.clusters = c("DAM", "Mono/Mac"),
  stage.lfc = 0.5, trem2.lfc = 0.58, padj = 0.05,
  stage.alpha = 0.1, # results() default used for the original Late vs Early tables
  trem2.min.cells = 20L, trem2.min.pairs = 2L,
  trem2.min.count = 10L, trem2.min.samples = 2L,
  go.min.trem2 = 10L, go.ont = "ALL", go.qvalue = 0.2,
  n.label = 5L, seed = 123L,
  make.markers = TRUE, run.go = TRUE
)

celltype.order = c("Homeostatic/early activated microglia", "DAM", "BAM",
                   "Mono/Mac", "Neutrophil", "MHC-II+ DC", "mregDC")
celltype.colors = c("Homeostatic/early activated microglia" = "#4E79A7",
                    "DAM" = "#E15759", "BAM" = "#59A14F", "Mono/Mac" = "#F28E2B",
                    "Neutrophil" = "#D4A72C", "MHC-II+ DC" = "#76B7B2",
                    "mregDC" = "#D45087", "Tumor" = "#969696", "NK" = "#2F2F2F")
micro.order = c("Homeostatic", "pre-DAM", "DAM1", "DAM2")
micro.colors = setNames(c("#4E79A7", "#59A14F", "#F28E2B", "#E15759"), micro.order)

myeloid.genes = c("P2ry12", "Tmem119", "Sall1", "Siglech", "Gpr34",
  "Cst7", "Apoc1", "Ccl12", "Bst2", "Ifi27l2a", "Lyve1", "Cd163", "Mrc1", "Folr2",
  "Cd209f", "Arg1", "Spp1", "Pf4", "Ms4a7", "Ccr2", "Ly6c2", "S100a8", "S100a9",
  "Retnlg", "Cxcr2", "Mmp9", "H2-Aa", "Cd74", "Flt3", "Itgax", "Clec10a", "Ccr7",
  "Fscn1", "Ccl22", "Cd274", "Pdcd1lg2")
micro.genes = c("P2ry12", "Tmem119", "Sall1", "Gpr34", "Mertk", "Atf3", "Egr1", "Cd83",
  "Tnf", "Ccl3", "Ifi27l2a", "Ifitm3", "Isg15", "Irf7", "Bst2", "Spp1", "Apoe",
  "Hmox1", "Fabp5", "Gpnmb")
expression.genes = c("Trem2", "Tfrc")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
samples = read.csv(text = "sample,stage
S1,Sham
S2,Sham
E1,Early
E2,Early
L1,Late
L2,Late
", stringsAsFactors = FALSE, na.strings = "")
stopifnot(!anyDuplicated(samples$sample), all(samples$stage %in% c("Sham", "Early", "Late")))
sample.stage = setNames(samples$stage, samples$sample)
set.seed(cfg$seed)
dir.create(cfg$output.dir, recursive = TRUE, showWarnings = FALSE)

safe_name = function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_table = function(x, file) write.csv(x, file, row.names = FALSE, na = "")
new_dir = function(name) {
  root = file.path(cfg$output.dir, name)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  # A separate directory per execution avoids stale results after a skipped test.
  out = tempfile(paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = root)
  dir.create(out)
  out
}
finish_run = function(out) {
  saveRDS(cfg, file.path(out, "settings.rds"))
  save_table(samples, file.path(out, "sample_metadata.csv"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  message("Results: ", normalizePath(out))
}
cell_ids = function(x) {
  ids = if ("cell" %in% names(x[[]])) as.character(x$cell) else colnames(x)
  if (anyNA(ids) || any(ids == "") || anyDuplicated(ids)) stop("Missing or duplicate cell identifiers.")
  ids
}
add_stage = function(x) {
  if (!"orig.ident" %in% names(x[[]])) stop("orig.ident is required.")
  unknown = setdiff(unique(as.character(x$orig.ident)), samples$sample)
  if (length(unknown)) stop("Samples absent from sample_metadata.csv: ", paste(unknown, collapse = ", "))
  x$stage = factor(unname(sample.stage[as.character(x$orig.ident)]), levels = c("Sham", "Early", "Late"))
  x
}
# Use original objects in the R session, or read their saved RDS files.
input_object = function(key) {
  object.names = c(all = "clean.seurat", myeloid = "myeloid.seurat.clean",
                   microglia = "microglia.seurat.clean")
  name = object.names[[key]]
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    x = get(name, envir = .GlobalEnv)
  } else {
    path = file.path(cfg$input.dir, cfg$inputs[[key]])
    if (!file.exists(path)) stop("Load the original object ", name, " or provide ", path)
    x = readRDS(path)
  }
  if (!inherits(x, "Seurat")) stop("Input is not a Seurat object: ", name)
  add_stage(x)
}
read_object = function(name) {
  p = file.path(cfg$output.dir, "objects", paste0(name, ".rds"))
  if (!file.exists(p)) stop("Run 01_annotation.R first: ", p)
  readRDS(p)
}
raw_counts = function(x) {
  if (!"RNA" %in% SeuratObject::Assays(x)) stop("RNA assay is required.")
  assay = x[["RNA"]]
  layers = grep("^counts($|\\.)", SeuratObject::Layers(assay, search = NA), value = TRUE)
  if (!length(layers)) stop("Raw RNA counts are required.")
  if (length(layers) > 1) {
    assay = SeuratObject::JoinLayers(assay, layers = "counts", new = "counts")
    layers = "counts"
  }
  m = SeuratObject::LayerData(assay, layer = layers[1])
  if (!all(colnames(x) %in% colnames(m))) stop("RNA counts are missing cells.")
  m[, colnames(x), drop = FALSE]
}
pseudobulk = function(m, groups, levels) {
  stopifnot(length(groups) == ncol(m), !anyNA(groups))
  n = tabulate(match(groups, levels), nbins = length(levels))
  if (any(n == 0)) stop("Pseudobulk group has no cells: ", paste(levels[n == 0], collapse = ", "))
  design = Matrix::sparseMatrix(i = seq_along(groups), j = match(groups, levels),
                                x = 1, dims = c(length(groups), length(levels)))
  pb = as.matrix(m %*% design)
  colnames(pb) = levels
  if (anyNA(pb) || any(!is.finite(pb)) || any(pb < 0) ||
      any(abs(pb - round(pb)) > 1e-8) || any(pb > .Machine$integer.max)) stop("Invalid raw pseudobulk counts.")
  storage.mode(pb) = "integer"
  pb
}
format_de = function(res, cluster, comparison, cutoff, stage = NA_character_) {
  tab = as.data.frame(res)
  tab$gene = rownames(tab)
  valid = !is.na(tab$padj) & is.finite(tab$log2FoldChange)
  tab$sig = ifelse(valid, "non-sig", "not tested/filtered")
  tab$sig[valid & tab$padj <= cfg$padj & tab$log2FoldChange >= cutoff] = "up"
  tab$sig[valid & tab$padj <= cfg$padj & tab$log2FoldChange <= -cutoff] = "down"
  tab$minus_log10_padj = -log10(pmax(tab$padj, 1e-300))
  tab$cluster = cluster; tab$comparison = comparison; tab$stage = stage
  tab = tab[order(tab$padj, -abs(tab$log2FoldChange), na.last = TRUE), ]
  rownames(tab) = NULL
  tab[, c("comparison", "cluster", "stage", "gene", setdiff(names(tab), c("comparison", "cluster", "stage", "gene")))]
}
volcano = function(tab, cutoff, title, xlab) {
  tab = tab[!is.na(tab$padj) & is.finite(tab$log2FoldChange), ]
  labels = tab %>% filter(sig %in% c("up", "down")) %>%
    group_by(sig) %>% slice_min(padj, n = cfg$n.label, with_ties = FALSE) %>% ungroup()
  ggplot(tab, aes(log2FoldChange, minus_log10_padj, color = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    ggrepel::geom_text_repel(data = labels, aes(label = gene), size = 3, max.overlaps = Inf,
                           seed = cfg$seed, min.segment.length = 0, show.legend = FALSE) +
    geom_vline(xintercept = c(-cutoff, cutoff), linetype = "dashed", color = "gray55") +
    geom_hline(yintercept = -log10(cfg$padj), linetype = "dashed", color = "gray55") +
    scale_color_manual(values = c(up = "#E15759", down = "#4E79A7", `non-sig` = "gray80")) +
    theme_classic() + theme(legend.position = "none") +
    labs(title = title, x = xlab, y = "-log10 adjusted P value")
}
save_de = function(fit, tab, pb, coldata, out, stem) {
  save_table(tab, file.path(out, paste0("DESeq2_", stem, ".csv")))
  save_table(tab[tab$sig %in% c("up", "down"), ], file.path(out, paste0("DEGs_", stem, ".csv")))
  write.csv(pb, file.path(out, paste0("Pseudobulk_counts_", stem, ".csv")))
  coldata$size_factor = DESeq2::sizeFactors(fit$dds)
  write.csv(coldata, file.path(out, paste0("Samples_", stem, ".csv")))
  write.csv(DESeq2::counts(fit$dds, normalized = TRUE), file.path(out, paste0("Normalized_counts_", stem, ".csv")))
  saveRDS(fit, file.path(out, paste0("DESeq2_model_", stem, ".rds")))
}
go_analysis = function(tab, out, stem, background = NULL, minimum = 1L) {
  tables = list(); summaries = list()
  for (direction in c("up", "down")) {
    genes = unique(tab$gene[tab$sig == direction])
    summary = data.frame(cluster = tab$cluster[1], stage = tab$stage[1], comparison = tab$comparison[1],
      direction = direction, n_DEGs = length(genes), n_terms = 0L,
      background = if (is.null(background)) "Default GO-annotated genes" else "Non-missing DESeq2 padj",
      status = "skipped: too few DEGs")
    if (!cfg$run.go) summary$status = "disabled in config"
    if (cfg$run.go && length(genes) >= minimum) {
      args = list(gene = genes, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL",
                  ont = cfg$go.ont, pvalueCutoff = cfg$padj, pAdjustMethod = "BH",
                  qvalueCutoff = cfg$go.qvalue)
      if (!is.null(background)) args$universe = background
      ego = tryCatch(do.call(clusterProfiler::enrichGO, args), error = function(e) e)
      if (inherits(ego, "error")) {
        summary$status = paste("failed:", conditionMessage(ego))
      } else {
        gt = if (is.null(ego)) data.frame() else as.data.frame(ego)
        summary$n_terms = nrow(gt)
        summary$status = if (nrow(gt)) "tested" else "no enriched terms"
        if (nrow(gt)) {
          gt$cluster = tab$cluster[1]; gt$stage = tab$stage[1]
          gt$comparison = tab$comparison[1]; gt$direction = direction
          tables[[direction]] = gt
          save_table(gt, file.path(out, paste0("GO_", stem, "_", direction, ".csv")))
          bp = if ("ONTOLOGY" %in% names(gt)) gt[gt$ONTOLOGY == "BP", ] else gt
          bp = head(bp[order(bp$p.adjust), ], 8)
          if (nrow(bp)) {
            bp$Description = factor(bp$Description, levels = rev(unique(bp$Description)))
            p = ggplot(bp, aes(-log10(pmax(p.adjust, 1e-300)), Description)) +
              geom_col(fill = if (direction == "up") "#E15759" else "#4E79A7") +
              scale_y_discrete(labels = function(x) stringr::str_wrap(x, 45)) + theme_classic() +
              labs(title = paste(stem, direction), x = "-log10 adjusted P value", y = NULL)
            ggsave(file.path(out, paste0("GO_", stem, "_", direction, ".pdf")), p,
                   width = 7, height = 2 + 0.35 * nrow(bp))
          }
        }
      }
    }
    summaries[[direction]] = summary
  }
  list(table = bind_rows(tables), summary = bind_rows(summaries))
}
