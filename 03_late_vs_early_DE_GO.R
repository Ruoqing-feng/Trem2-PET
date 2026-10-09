source("00_settings.R")
suppressPackageStartupMessages(library(DESeq2))
out = new_dir("Late_vs_Early")
go.dir = file.path(out, "GO")
dir.create(go.dir)
obj = read_object("myeloid")
sample.order = samples$sample[samples$stage %in% c("Early", "Late")]
stopifnot(sum(sample.stage[sample.order] == "Early") >= 2, sum(sample.stage[sample.order] == "Late") >= 2)
de.results = list(); de.summary = list(); go.results = list(); go.summary = list(); plots = list()

for (cl in cfg$stage.clusters) {
  message("Late vs Early: ", cl)
  cells = colnames(obj)[obj$sub4 %in% cl & obj$orig.ident %in% sample.order]
  n = table(factor(as.character(obj$orig.ident[match(cells, colnames(obj))]), levels = sample.order))
  stem = safe_name(cl)
  coldata = data.frame(sample = sample.order,
    condition = factor(unname(sample.stage[sample.order]), levels = c("Early", "Late")),
    n_cells = as.integer(n), row.names = sample.order)
  save_table(coldata, file.path(out, paste0("Cell_counts_", stem, ".csv")))
  summary = data.frame(cluster = cl, n_samples = sum(n > 0), n_up = NA_integer_, n_down = NA_integer_, status = "not run")
  if (any(n == 0)) {
    summary$status = paste("skipped: no cells in", paste(names(n)[n == 0], collapse = ", "))
    de.summary[[cl]] = summary
    next
  }
  sub = subset(obj, cells = cells)
  counts = raw_counts(sub)
  # Same raw sums as AggregateExpression(group.by = 'orig.ident').
  pb = pseudobulk(counts, as.character(sub$orig.ident), sample.order)
  fit = tryCatch({
    dds = DESeqDataSetFromMatrix(pb, coldata, design = ~condition)
    dds = DESeq(dds)
    # No extra gene prefilter. Keep the original default independent-filtering alpha.
    res = results(dds, contrast = c("condition", "Late", "Early"), alpha = cfg$stage.alpha)
    list(dds = dds, res = res)
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    summary$status = paste("failed:", conditionMessage(fit))
    de.summary[[cl]] = summary
    next
  }
  tab = format_de(fit$res, cl, "Late_vs_Early", cfg$stage.lfc)
  de.results[[cl]] = tab
  save_de(fit, tab, pb, coldata, out, stem)
  summary$n_up = sum(tab$sig == "up"); summary$n_down = sum(tab$sig == "down")
  summary$status = "tested"; de.summary[[cl]] = summary
  p = volcano(tab, cfg$stage.lfc, cl, "log2 fold change (Late stage vs early stage)")
  plots[[cl]] = p
  ggsave(file.path(out, paste0("Volcano_", stem, ".pdf")), p, width = 5, height = 5)
  go = go_analysis(tab, go.dir, stem) # Original annotated-gene background; no custom universe.
  go.results[[cl]] = go$table; go.summary[[cl]] = go$summary
}
save_table(bind_rows(de.summary), file.path(out, "DE_summary.csv"))
saveRDS(de.results, file.path(out, "DESeq2_results.rds"))
all.de = bind_rows(de.results)
if (nrow(all.de)) save_table(all.de[all.de$sig %in% c("up", "down"), ], file.path(out, "DEGs_all_comparisons.csv"))
all.go = bind_rows(go.results)
saveRDS(go.results, file.path(go.dir, "GO_results.rds"))
save_table(bind_rows(go.summary), file.path(go.dir, "GO_summary.csv"))
if (nrow(all.go)) save_table(all.go, file.path(go.dir, "GO_all_comparisons.csv"))
if (length(plots)) ggsave(file.path(out, "Volcano_Late_vs_Early_panel.pdf"), wrap_plots(plots, ncol = 2),
                        width = 10, height = 5 * ceiling(length(plots) / 2))

# Named terms reproduce the selections in new_analysis.R and Figure1_GO_barplots.R.
# Report missing/non-significant selections; never substitute a term from another comparison.
if (cfg$run.go) {
  picks = read.csv(text = "cluster,direction,ID,Description
DAM,up,,leukocyte chemotaxis
DAM,up,,response to lipopolysaccharide
DAM,up,,positive regulation of innate immune response
DAM,up,,tumor necrosis factor production
DAM,up,,pattern recognition receptor signaling pathway
DAM,up,,canonical NF-kappaB signal transduction
DAM,up,,chemokine-mediated signaling pathway
DAM,down,,negative regulation of inflammasome-mediated signaling pathway
Mono/Mac,up,,response to transforming growth factor beta
Mono/Mac,up,,negative regulation of inflammatory response
Mono/Mac,up,,negative regulation of T cell proliferation
Mono/Mac,up,,arginine metabolic process
Mono/Mac,down,,canonical NF-kappaB signal transduction
Mono/Mac,down,,tumor necrosis factor production
Mono/Mac,down,,inflammasome-mediated signaling pathway
Mono/Mac,down,,antigen processing and presentation
MHC-II+ DC,up,,leukocyte chemotaxis
MHC-II+ DC,up,,myeloid leukocyte differentiation
MHC-II+ DC,up,,phagocytosis
MHC-II+ DC,up,,response to transforming growth factor beta
MHC-II+ DC,up,,negative regulation of cytokine production
MHC-II+ DC,down,,response to interferon-beta
MHC-II+ DC,down,,response to type II interferon
MHC-II+ DC,down,,antigen processing and presentation
MHC-II+ DC,down,,positive regulation of leukocyte cell-cell adhesion
Homeostatic/early activated microglia,up,,leukocyte chemotaxis
Homeostatic/early activated microglia,up,,response to lipopolysaccharide
Homeostatic/early activated microglia,up,,response to oxidative stress
Homeostatic/early activated microglia,up,,positive regulation of angiogenesis
Homeostatic/early activated microglia,up,,antigen processing and presentation
Homeostatic/early activated microglia,down,,regulation of synapse organization
Homeostatic/early activated microglia,down,,glial cell migration
", stringsAsFactors = FALSE, na.strings = "")
  picks$found = FALSE
  selected = list()
  for (i in seq_len(nrow(picks))) {
    if (!nrow(all.go)) next
    hit = all.go$cluster == picks$cluster[i] & all.go$direction == picks$direction[i] & all.go$ONTOLOGY == "BP"
    hit = hit & if (!is.na(picks$ID[i])) all.go$ID == picks$ID[i] else all.go$Description == picks$Description[i]
    picks$found[i] = any(hit)
    if (any(hit)) selected[[as.character(i)]] = all.go[hit, ]
  }
  save_table(picks, file.path(go.dir, "Selected_terms_check.csv"))
  selected = bind_rows(selected)
  if (nrow(selected)) {
    selected = unique(selected)
    save_table(selected, file.path(go.dir, "GO_selected_terms.csv"))
    for (cl in unique(selected$cluster)) {
      t = selected[selected$cluster == cl, ]
      t$value = -log10(pmax(t$p.adjust, 1e-300)) * ifelse(t$direction == "up", 1, -1)
      t$Description = factor(t$Description, levels = unique(t$Description[order(t$value)]))
      p = ggplot(t, aes(value, Description, fill = direction)) + geom_col() +
        geom_vline(xintercept = 0, color = "grey60") +
        scale_fill_manual(values = c(up = "#E15759", down = "#4E79A7"),
          breaks = c("up", "down"), labels = c("Higher in late stage", "Higher in early stage")) +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, 45)) + theme_classic() +
        labs(title = cl, x = "Signed -log10 adjusted P value", y = NULL, fill = NULL)
      ggsave(file.path(go.dir, paste0("GO_selected_", safe_name(cl), ".pdf")), p,
             width = 8, height = 2 + 0.4 * nrow(t))
    }
  }
}
print(bind_rows(de.summary))
finish_run(out)
