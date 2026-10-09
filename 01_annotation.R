# Inputs must retain their original labels, especially microglia.seurat.clean$subtype.
source("00_settings.R")
out = file.path(cfg$output.dir, "annotation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$output.dir, "objects"), recursive = TRUE, showWarnings = FALSE)

# Annotate by stable cell identifiers. Do not infer labels from a new clustering run.
myeloid = input_object("myeloid")
microglia = input_object("microglia")
stopifnot("sub2" %in% names(myeloid[[]]), "subtype" %in% names(microglia[[]]))
map = read.csv(text = "source,source_label,final_label,keep,reason
myeloid_sub2,TAM,BAM,TRUE,
myeloid_sub2,BAM,BAM,TRUE,
myeloid_sub2,Monocyte,mregDC,TRUE,
myeloid_sub2,mregDC,mregDC,TRUE,
myeloid_sub2,Neutrophil,Neutrophil,TRUE,
microglia_subtype,Homeostatic,Homeostatic,TRUE,
microglia_subtype,pre-DAM,Early activated microglia,TRUE,
microglia_subtype,DAM1,DAM1,TRUE,
microglia_subtype,DAM3,DAM2,TRUE,
microglia_subtype,DAM2,Arg1+ Spp1+ macrophages,TRUE,
microglia_subtype,DAM4,MHC-II+ DC,TRUE,
microglia_subtype,DAM5,Ccr2+ inflammatory monocytes,TRUE,
microglia_subtype,DAM6,,FALSE,Putative doublets with strong tumor and myeloid marker co-expression
", stringsAsFactors = FALSE, na.strings = "")
stopifnot(!anyDuplicated(paste(map$source, map$source_label)))
md = myeloid[[]]
md$cell = cell_ids(myeloid)
md$source_label = as.character(md$sub2)
md$source = "myeloid_sub2"
mi = match(md$cell, cell_ids(microglia))
has.micro = !is.na(mi)
md$source_label[has.micro] = as.character(microglia$subtype[mi[has.micro]])
md$source[has.micro] = "microglia_subtype"
ix = match(paste(md$source, md$source_label), paste(map$source, map$source_label))
if (anyNA(ix)) stop("Unmapped source annotations: ",
  paste(unique(paste(md$source[is.na(ix)], md$source_label[is.na(ix)])), collapse = ", "),
  ". Check that input objects precede final annotation.")
md$keep = map$keep[ix]
md$reason = map$reason[ix]
md$celltype_fine = map$final_label[ix]
myeloid = subset(myeloid, cells = rownames(md)[md$keep])
myeloid$cell = md$cell[md$keep]
myeloid$sub2 = md$celltype_fine[md$keep]
myeloid$sub3 = as.character(myeloid$sub2)
myeloid$sub3[myeloid$sub3 %in% c("Homeostatic", "Early activated microglia")] = "Homeostatic/early activated microglia"
myeloid$sub3[myeloid$sub3 %in% c("DAM1", "DAM2")] = "DAM"
myeloid$sub4 = as.character(myeloid$sub3)
myeloid$sub4[myeloid$sub4 %in% c("Arg1+ Spp1+ macrophages", "Ccr2+ inflammatory monocytes")] = "Mono/Mac"
stopifnot(all(myeloid$sub4 %in% celltype.order))
Idents(myeloid) = "sub4"

micro.cells = colnames(microglia)[cell_ids(microglia) %in%
  myeloid$cell[myeloid$sub2 %in% c("Homeostatic", "Early activated microglia", "DAM1", "DAM2")]]
if (!length(micro.cells)) stop("No retained microglial cells.")
microglia = subset(microglia, cells = micro.cells)
microglia$cell = cell_ids(microglia)
microglia$subtype = as.character(myeloid$sub2[match(microglia$cell, myeloid$cell)])
microglia$subtype[microglia$subtype == "Early activated microglia"] = "pre-DAM"
Idents(microglia) = factor(microglia$subtype, levels = micro.order)

all.cells = input_object("all")
stopifnot("major" %in% names(all.cells[[]]))
all.cells$cell = cell_ids(all.cells)
if (!all(myeloid$cell %in% all.cells$cell)) stop("Myeloid barcodes are absent from the all-cell input.")
all.md = all.cells[[]]
myeloid.match = match(all.md$cell, myeloid$cell)
is.myeloid = !is.na(myeloid.match)
is.other = all.md$major %in% c("Tumor", "NK")
if (any(is.myeloid & is.other)) stop("Conflicting major and myeloid identities.")
unknown.major = !all.md$major %in% c("Mac/microglia", "Monocyte", "Myeloid", "Tumor", "NK")
if (any(unknown.major)) stop("Unrecognized all-cell major annotations: ", paste(unique(all.md$major[unknown.major]), collapse = ", "))
all.md$keep = is.myeloid | is.other
all.md$reason = NA_character_
all.md$reason[!all.md$keep] = "Not retained in the supplied curated myeloid subset"
excluded.match = match(all.md$cell, md$cell[!md$keep])
ii = which(!is.na(excluded.match))
all.md$reason[ii] = md$reason[!md$keep][excluded.match[ii]]
all.md$celltype = as.character(all.md$major)
all.md$celltype_fine = as.character(all.md$major)
all.md$celltype[is.myeloid] = as.character(myeloid$sub4[myeloid.match[is.myeloid]])
all.md$celltype_fine[is.myeloid] = as.character(myeloid$sub2[myeloid.match[is.myeloid]])
all.md$celltype[!all.md$keep] = NA_character_
all.md$celltype_fine[!all.md$keep] = NA_character_
all.cells = subset(all.cells, cells = rownames(all.md)[all.md$keep])
all.cells$myeloid_subtype2 = all.md$celltype[all.md$keep]
all.cells$celltype_fine = all.md$celltype_fine[all.md$keep]
all.cells$major = ifelse(all.cells$myeloid_subtype2 %in% celltype.order, "Myeloid", all.cells$myeloid_subtype2)

# Reuse saved coordinates; never rerun normalization, PCA, graph clustering or UMAP.
use_coordinates = function(x, file, object.name) {
  path = file.path(cfg$input.dir, file)
  coords = NULL
  coordinate.source = "Inherited input coordinates; no UMAP recomputation"
  if (exists(object.name, envir = .GlobalEnv, inherits = FALSE)) {
    original = get(object.name, envir = .GlobalEnv)
    coords = Embeddings(original, "umap")[colnames(original), , drop = FALSE]
    rownames(coords) = cell_ids(original)
    coordinate.source = object.name
  } else if (file.exists(path)) {
    coords = readRDS(path)
    coordinate.source = file
  }
  if (!is.null(coords)) {
    if (!is.matrix(coords) || anyDuplicated(rownames(coords)) ||
        !all(cell_ids(x) %in% rownames(coords))) stop("Invalid or incomplete embedding: ", path)
    coords = coords[cell_ids(x), , drop = FALSE]
    rownames(coords) = colnames(x)
    colnames(coords) = paste0("UMAP_", seq_len(ncol(coords)))
    x[["umap"]] = CreateDimReducObject(coords, key = "UMAP_", assay = DefaultAssay(x))
    attr(x, "coordinate_source") = coordinate.source
  } else {
    if (!"umap" %in% SeuratObject::Reductions(x)) stop("No saved UMAP in input object.")
    attr(x, "coordinate_source") = "Inherited input coordinates; no UMAP recomputation"
  }
  x
}
all.cells = use_coordinates(all.cells, cfg$all.umap, "all.seurat.clean")
microglia = use_coordinates(microglia, cfg$microglia.umap, "real.microglia.seurat.clean")
save_table(data.frame(object = c("all", "microglia"),
  coordinate_source = c(attr(all.cells, "coordinate_source"), attr(microglia, "coordinate_source"))),
  file.path(out, "embedding_sources.csv"))

save_table(md[, c("cell", "orig.ident", "stage", "source", "source_label", "celltype_fine", "keep", "reason")],
           file.path(out, "annotation_provenance.csv"))
ledger = all.md[, c("cell", "orig.ident", "stage", "major", "celltype_fine", "celltype", "keep", "reason")]
save_table(ledger, file.path(out, "cell_annotations.csv"))
save_table(ledger[!ledger$keep, ], file.path(out, "excluded_cells.csv"))
save_table(as.data.frame(table(myeloid$orig.ident, myeloid$sub4)), file.path(out, "cell_counts_by_sample.csv"))
save_table(map, file.path(out, "annotation_mapping.csv"))
for (name in c("all.cells", "myeloid", "microglia")) {
  saveRDS(get(name), file.path(cfg$output.dir, "objects", paste0(name, ".rds")))
}
finish_run(out)
