source("00_settings.R")
out = new_dir("gene_detection")
obj = read_object("all.cells")
counts = raw_counts(obj)
if (!all(expression.genes %in% rownames(counts))) stop("Missing RNA gene(s): ", paste(setdiff(expression.genes, rownames(counts)), collapse = ", "))
md = obj[[]]
populations = unique(c(celltype.order, "Tumor", "NK", "Myeloid"))
rows = list()
for (gene in expression.genes) {
  v = as.numeric(counts[gene, ])
  if (anyNA(v) || any(!is.finite(v)) || any(v < 0) || any(abs(v - round(v)) > 1e-8)) stop("Invalid raw RNA counts: ", gene)
  for (pop in populations) {
    in.pop = if (pop == "Myeloid") md$major %in% "Myeloid" else md$myeloid_subtype2 %in% pop
    for (s in samples$sample) {
      ii = in.pop & md$orig.ident %in% s
      rows[[paste(gene, pop, s)]] = data.frame(gene = gene, population = pop, sample = s,
        stage = unname(sample.stage[s]), n_cells = sum(ii), n_positive = sum(v[ii] > 0),
        percent_positive = if (any(ii)) 100 * mean(v[ii] > 0) else NA_real_)
    }
  }
}
by.sample = bind_rows(rows)
by.sample$stage = factor(by.sample$stage, levels = c("Sham", "Early", "Late"))
by.stage = by.sample %>% group_by(gene, population, stage) %>%
  summarise(n_samples = sum(n_cells > 0), n_cells = sum(n_cells), n_positive = sum(n_positive),
    percent_pooled = ifelse(n_cells > 0, 100 * n_positive / n_cells, NA_real_),
    percent_sample_mean = if (all(is.na(percent_positive))) NA_real_ else mean(percent_positive, na.rm = TRUE),
    percent_sample_sd = sd(percent_positive, na.rm = TRUE), .groups = "drop")
save_table(by.sample, file.path(out, "Detection_by_sample.csv"))
save_table(by.stage, file.path(out, "Detection_by_stage.csv"))
for (gene in expression.genes) {
  plot.tab = by.stage %>% filter(.data$gene == .env$gene, population != "Myeloid")
  plot.tab$population = factor(plot.tab$population, levels = c(celltype.order, "Tumor", "NK"))
  p = ggplot(plot.tab, aes(population, percent_pooled / 100, fill = stage)) +
    geom_col(position = "dodge", na.rm = TRUE) +
    scale_fill_manual(values = c(Sham = "#7F7F7F", Early = "#4E79A7", Late = "#E15759")) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = gene, x = NULL, y = "Cells with detectable RNA", fill = NULL)
  ggsave(file.path(out, paste0(gene, "_detection_by_stage.pdf")), p, width = 11, height = 5)
}
writeLines(c("Detection = raw RNA count >0. Absent populations have NA percentage, not zero.",
  "percent_pooled weights cells; percent_sample_mean weights animals equally (only samples with cells).",
  "Myeloid is an additional aggregate of its subtypes and must not be added to subtype counts.",
  "Inspect n_cells before interpreting percentages, particularly rare tumor-labelled cells in sham."),
  file.path(out, "percentage_definitions.txt"))
print(by.stage, n = Inf)
finish_run(out)
