source("00_settings.R")
suppressPackageStartupMessages(library(DESeq2))
out = new_dir("Trem2_positive_vs_negative")
go.dir = file.path(out, "GO")
dir.create(go.dir)
obj = read_object("myeloid")
de.results = list(); de.summary = list(); qc = list(); go.results = list(); go.summary = list(); plots = list()

for (cl in cfg$trem2.clusters) {
  for (st in c("Early", "Late")) {
    key = paste(cl, st, sep = " | "); stem = paste(safe_name(cl), st, sep = "_")
    sample.order = samples$sample[samples$stage == st]
    cells = colnames(obj)[obj$sub4 %in% cl & obj$orig.ident %in% sample.order]
    info = expand.grid(status = c("negative", "positive"), sample = sample.order, stringsAsFactors = FALSE)
    info$cluster = cl; info$stage = st; info$n_cells = 0L
    info$median_RNA_counts = NA_real_; info$median_RNA_features = NA_real_
    if (length(cells)) {
      sub = subset(obj, cells = cells)
      m = raw_counts(sub)
      if (!"Trem2" %in% rownames(m)) stop("Trem2 absent from RNA counts.")
      v = as.numeric(m["Trem2", ])
      if (anyNA(v) || any(v < 0) || any(!is.finite(v)) || any(abs(v - round(v)) > 1e-8)) stop("Invalid Trem2 raw counts.")
      md = sub[[]]
      md$status = ifelse(v > 0, "positive", "negative")
      md$RNA_counts = Matrix::colSums(m); md$RNA_features = Matrix::colSums(m > 0)
      for (j in seq_len(nrow(info))) {
        ii = md$orig.ident %in% info$sample[j] & md$status %in% info$status[j]
        info$n_cells[j] = sum(ii)
        if (any(ii)) {
          info$median_RNA_counts[j] = median(md$RNA_counts[ii])
          info$median_RNA_features[j] = median(md$RNA_features[ii])
        }
      }
    }
    total = ave(info$n_cells, info$sample, FUN = sum)
    info$percentage = ifelse(total > 0, 100 * info$n_cells / total, NA_real_)
    paired = sample.order[vapply(sample.order, function(s) all(info$n_cells[info$sample == s] >= cfg$trem2.min.cells), logical(1))]
    info$included = info$sample %in% paired
    qc[[key]] = info
    summary = data.frame(cluster = cl, stage = st, n_pairs = length(paired),
      samples_used = paste(paired, collapse = ", "), excluded_samples = paste(setdiff(sample.order, paired), collapse = ", "),
      n_genes_test_input = NA_integer_, n_up = NA_integer_, n_down = NA_integer_, status = "not run")
    if (length(paired) < cfg$trem2.min.pairs) {
      summary$status = "skipped: insufficient complete animal pairs"
      de.summary[[key]] = summary
      next
    }
    cd = info[info$included, c("sample", "status", "n_cells")]
    cd$sample = factor(cd$sample, levels = paired)
    cd$status = factor(cd$status, levels = c("negative", "positive"))
    rownames(cd) = paste(cd$sample, cd$status, sep = "_")
    use = md$orig.ident %in% paired
    pb = pseudobulk(m[, use, drop = FALSE], paste(md$orig.ident[use], md$status[use], sep = "_"), rownames(cd))
    keep = rownames(pb) != "Trem2" & rowSums(pb >= cfg$trem2.min.count) >= cfg$trem2.min.samples
    save_table(data.frame(gene = rownames(pb), included = keep), file.path(out, paste0("Gene_filter_", stem, ".csv")))
    write.csv(pb, file.path(out, paste0("Pseudobulk_counts_", stem, ".csv")))
    if (!any(keep) || any(colSums(pb[keep, , drop = FALSE]) == 0)) {
      summary$status = "skipped: no usable counts after filtering"
      de.summary[[key]] = summary
      next
    }
    fit = tryCatch({
      dds = DESeqDataSetFromMatrix(pb[keep, , drop = FALSE], cd, design = ~sample + status)
      dds = DESeq(dds)
      res = results(dds, contrast = c("status", "positive", "negative"), alpha = cfg$padj)
      list(dds = dds, res = res)
    }, error = function(e) e)
    if (inherits(fit, "error")) {
      summary$status = paste("failed:", conditionMessage(fit)); de.summary[[key]] = summary
      next
    }
    tab = format_de(fit$res, cl, "Trem2_positive_vs_negative", cfg$trem2.lfc, st)
    de.results[[key]] = tab
    summary$n_genes_test_input = nrow(tab)
    summary$n_up = sum(tab$sig == "up"); summary$n_down = sum(tab$sig == "down")
    summary$status = "tested"; de.summary[[key]] = summary
    save_de(fit, tab, pb, cd, out, stem)
    p = volcano(tab, cfg$trem2.lfc, key, "log2 fold change (Trem2 detected vs undetected)")
    plots[[key]] = p
    ggsave(file.path(out, paste0("Volcano_", stem, ".pdf")), p, width = 5.5, height = 5)
    universe = tab$gene[!is.na(tab$padj)]
    save_table(data.frame(gene = universe), file.path(go.dir, paste0("Universe_", stem, ".csv")))
    go = go_analysis(tab, go.dir, stem, background = universe, minimum = cfg$go.min.trem2)
    go.results[[key]] = go$table; go.summary[[key]] = go$summary
  }
}
save_table(bind_rows(de.summary), file.path(out, "DE_summary.csv"))
save_table(bind_rows(qc), file.path(out, "sample_QC.csv"))
saveRDS(de.results, file.path(out, "DESeq2_results.rds"))
all.de = bind_rows(de.results)
if (nrow(all.de)) save_table(all.de[all.de$sig %in% c("up", "down"), ], file.path(out, "DEGs_all_comparisons.csv"))
saveRDS(go.results, file.path(go.dir, "GO_results.rds"))
save_table(bind_rows(go.summary), file.path(go.dir, "GO_summary.csv"))
all.go = bind_rows(go.results)
if (nrow(all.go)) save_table(all.go, file.path(go.dir, "GO_all_comparisons.csv"))
if (length(plots)) ggsave(file.path(out, "Volcano_Trem2_status_panel.pdf"), wrap_plots(plots, ncol = 2),
                        width = 11, height = 5 * ceiling(length(plots) / 2))
writeLines(c("Detected = raw RNA Trem2 counts >0; undetected = 0, not proven absence of protein.",
  "Associations can reflect detection depth and within-population composition; check sample_QC.csv.",
  "Two animals per stage: exploratory comparisons, not a stage-by-status interaction test.",
  "Trem2 is excluded from DESeq2 normalization, testing and GO.",
  "The log2FC cutoff is applied after the zero-effect Wald test."), file.path(out, "analysis_notes.txt"))
print(bind_rows(de.summary))
finish_run(out)
