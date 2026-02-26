###############################################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(fgsea)
  library(data.table)
  library(scCustomize)
})

set.seed(1)

# ------------------------------- paths ---------------------------------------
out_dir <- "G:/PhD_final/tables"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------- load data -----------------------------------
load("G:/PhD_final/result_obj_new.RData")                  # result_obj (Seurat)
load("G:/PhD_final/tRNA_miRNA_selected_raw_counts.RData")  # tRNA_miRNA_selected_raw_counts
load("G:/PhD_final/all_stringtie_selected.RData")          # all_stringtie_selected
load("D:/scRNA-seq/tRF_motif/cds_tx_seq.RData")            # cds_tx_seq
load("D:/scRNA-seq/tRF_motif/tx2gene.RData")               # tx2gene
load("D:/scRNA-seq/tRF_motif/utr_tx_seq.RData")            # utr_tx_seq
#load("D:/scRNA-seq/AZ_final_obj/filtered_DEG.RData")       # filtered_DE (optional)
load("D:/scRNA-seq/AZ_final_obj/filtered_DEG_abr.RData") #filtered_DE
load("G:/PhD_final/final_bulk_DGE.RData")                  # optional
load("G:/PhD_final/all_stringtie_raw_counts_counts.RData") # optional

result_obj@meta.data
unique(result_obj$final_population)

Idents(result_obj) <- result_obj$cluster_key_final
#result_obj$orig.ident
#DimPlot(subset(result_obj, idents = "g24.2"),reduction = "umap.d33.nn100.md0.3")
# 1. Identify the cell barcodes for the specific cluster (e.g., cluster "3")
cells_to_highlight <- WhichCells(result_obj, idents = "g24.2")
cells_to_highlight <- WhichCells(result_obj, idents = c("g9.1","g9.2","g9.3","g9.4"))
cells_to_highlight <- WhichCells(result_obj, idents = c("g32.1","g35.1","g35.2","g35.3","g35.4"))
cells_to_highlight <- WhichCells(result_obj, idents = "g22.1.3")
cells_to_highlight <- WhichCells(result_obj, idents = c("g15.1","g15.2"))
cells_to_highlight <- WhichCells(result_obj, idents = "g15.2")
cells_to_highlight <- WhichCells(result_obj, idents = c("g7.9","g50.4","g52.4","g7.2"))
cells_to_highlight <- WhichCells(result_obj, idents = "g14.4")
cells_to_highlight <- WhichCells(result_obj, idents = c("g50.3","g52.3","g7.6"))


cells_to_highlight <- WhichCells(result_obj, 
                                 idents = c("g1.1.1","g1.1.2","g1.1.3","g1.1.4","g1.2","g1.3",
                                            "g12.1.4","g12.2","g26.2","g26.3","g30.3","g39.2"))
cells_to_highlight <- WhichCells(result_obj, idents = c("g18.1.1","g18.1.2","g18.2"))


cells_to_highlight <- WhichCells(result_obj, idents = "g29.2")


cells_to_highlight <- WhichCells(result_obj, idents = c("g49.1","g49.2"))
cells_to_highlight <- WhichCells(result_obj, idents = "g38.1")
cells_to_highlight <- WhichCells(result_obj, idents = c("g3.1.1.3","g3.3"))

cells_to_highlight <- WhichCells(result_obj, idents = c("g7.1.1","g7.1.2","g7.1.3","g7.1.4",
                                                        "g7.10","g7.3","g7.4"))
cells_to_highlight <- WhichCells(result_obj, idents = "g7.4")



cells_to_highlight <- WhichCells(result_obj, idents = "g50.4")
cells_to_highlight <- WhichCells(result_obj, idents = "g21.3")


# 2. Plot with highlighting
DimPlot(result_obj, 
        cells.highlight = cells_to_highlight, 
        cols.highlight = "red", 
        cols = "grey",
        reduction = "umap.d33.nn100.md0.3")









################################################################################
#Add abracada cells:
#gene_marker <- "SMESG000049716.1"#CTSL2
#gene_marker <- "SMESG000002628.1" #PGRN
#innexin 8 (inx8)	dd_Smed_v6_1234_0_1 SMESG000080705.1
#innexin 9 (inx9)	dd_Smed_v6_692_0_1 SMESG000080788.1
gene_marker <- c("SMESG000080705.1","SMESG000080788.1")#inx8
#SMESG000040772.1 | MMP1
#load("G:/PhD_final/result_obj_new.RData")                  # result_obj: Seurat obj
FeaturePlot(result_obj, features = "SMESG000056326.1", reduction = "umap.d33.nn100.md0.3")
FeaturePlot(result_obj, features = "SMESG000040772.1", reduction = "umap.d33.nn100.md0.3")
library(scCustomize)
?FeaturePlot_scCustom
result_obj$condition_correct
FeaturePlot_scCustom(seurat_object = result_obj,
                     split.by = "condition_correct",
                     features = "SMESG000056326.1",
                     reduction = "umap.d33.nn100.md0.3",
                     num_columns = 4
)
FeaturePlot_scCustom(seurat_object = result_obj,
                     split.by = "condition_correct",
                     features = "SMESG000036051.1",
                     reduction = "umap.d33.nn100.md0.3",
                     num_columns = 4
)

FeaturePlot_scCustom(seurat_object = result_obj,
                     split.by = "condition_correct",
                     features = "SMESG000015210.1",
                     reduction = "umap.d33.nn100.md0.3",
                     num_columns = 4
)

FeaturePlot_scCustom(seurat_object = result_obj,
                     split.by = "condition_correct",
                     features = "SMESG000030475.1",
                     reduction = "umap.d33.nn100.md0.3",
                     num_columns = 4
)

stopifnot("cluster_key_final" %in% colnames(result_obj@meta.data))
gene_marker <- c("SMESG000049716.1","SMESG000080705.1", "SMESG000080788.1")  # inx8, inx9

expr_mat <- NULL
used_assay <- NULL

for (assay in c("RNA", "SCT")) {
  if (assay %in% Seurat::Assays(result_obj)) {
    rn <- rownames(result_obj[[assay]])
    if (all(gene_marker %in% rn)) {
      expr_mat <- Seurat::GetAssayData(result_obj, assay = assay, layer = "data")[gene_marker, , drop = FALSE]
      used_assay <- assay
      break
    }
  }
}

if (is.null(expr_mat)) {
  stop("Neither RNA nor SCT contains *both* inx8 and inx9 features: ", paste(gene_marker, collapse = ", "))
}

used_assay
dim(expr_mat)
library(Matrix)
library(dplyr)

stopifnot("cluster_key_final" %in% colnames(result_obj@meta.data))

cells <- colnames(expr_mat)

# ensure cells align to meta.data (robust to any mismatch)
cells <- intersect(cells, rownames(result_obj@meta.data))
expr_mat2 <- expr_mat[, cells, drop = FALSE]
meta <- result_obj@meta.data[cells, , drop = FALSE]

# make a small dense matrix: cells × genes (74403 × 2) is fine
m <- as.matrix(t(expr_mat2))  # rows = cells, cols = genes

df <- data.frame(
  cell = cells,
  cluster_key_final = meta$cluster_key_final,
  m,
  check.names = FALSE
)

cluster_summ <- df |>
  dplyr::group_by(cluster_key_final) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    dplyr::across(
      .cols = all_of(colnames(m)),
      .fns  = list(mean = ~mean(.x, na.rm = TRUE),
                   pct  = ~mean(.x > 0, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(.data[[paste0(colnames(m)[1], "_mean")]]))

cluster_summ
as.data.frame(cluster_summ)
# Pull expression vector from RNA (preferred) or SCT as fallback
# expr_vec <- NULL
# used_assay <- NULL
# for (assay in c("RNA", "SCT")) {
#   if (assay %in% Seurat::Assays(result_obj) && gene_marker %in% rownames(result_obj[[assay]])) {
#     expr_vec <- Seurat::GetAssayData(result_obj, assay = assay, layer = "data")[gene_marker, ]
#     used_assay <- assay
#     break
#   }
# }
# if (is.null(expr_vec)) {
#   stop("Marker gene not found in RNA or SCT assay rownames: ", gene_marker)
# }
# message("Using assay: ", used_assay)

# df_marker <- data.frame(
#   cell = names(expr_mat),
#   cluster_key_final = result_obj@meta.data[names(expr_mat), "cluster_key_final"],
#   expr = as.numeric(expr_mat),
#   stringsAsFactors = FALSE
# )
# 
# df_marker <- df_marker[!is.na(df_marker$cluster_key_final), ]
# 
# cluster_summ <- df_marker |>
#   dplyr::group_by(cluster_key_final) |>
#   dplyr::summarise(
#     n_cells = dplyr::n(),
#     mean_expr = mean(expr, na.rm = TRUE),
#     pct_expr  = mean(expr > 0, na.rm = TRUE),
#     .groups = "drop"
#   ) |>
#   dplyr::arrange(dplyr::desc(mean_expr))
# cluster_summ
#print(head(cluster_summ, 15))
#cluster_summ[cluster_summ$pct_expr>0.8,]


################################################################################

load("G:/PhD_final/tables/targets_list_noCDS.RData")
targets_list_noCDS

load("G:/PhD_final/tables/targets_list_withCDS.RData")
targets_list_withCDS



# inputs
trf_id <- "GCATCGGTGGTTCAGTGGTAGAATGCTCGCCT dd_Smed_g4_GLY-CCC_tRNA_2"
genes_noCDS  <- targets_list_noCDS[[trf_id]]
genes_withCDS  <- targets_list_withCDS[[trf_id]]


load("D:/scRNA-seq/AZ_final_obj/filtered_DEG_abr_new.RData")
filtered_DEG_abr_new <- filtered_DE
load("D:/scRNA-seq/AZ_final_obj/filtered_DEG_abr.RData")
filtered_DEG_abr <- filtered_DE



filtered_DEG_abr[filtered_DEG_abr$gene %in% genes_noCDS & filtered_DEG_abr$time=="72h" &
                   filtered_DEG_abr$direction=="overexpressed",]
filtered_DEG_abr[filtered_DEG_abr$gene %in% genes_withCDS & filtered_DEG_abr$time=="72h" &
                   filtered_DEG_abr$direction=="overexpressed",]

filtered_DEG_abr_new[filtered_DEG_abr_new$gene %in% genes_noCDS & filtered_DEG_abr_new$time=="72h" &
                       filtered_DEG_abr_new$direction=="overexpressed",]
filtered_DEG_abr_new[filtered_DEG_abr_new$gene %in% genes_withCDS & filtered_DEG_abr_new$time=="72h" &
                       filtered_DEG_abr_new$direction=="overexpressed",]

################################################################################

# choose assay/slot (adjust if needed)
assay_use <- Seurat::DefaultAssay(result_obj)   # usually "RNA"
slot_use  <- "data"                            # log-normalized expression

# keep only genes present in object
genes_use <- intersect(genes, rownames(result_obj))
length(genes_use)

# average expression per final_population
avg_list <- Seurat::AverageExpression(
  object   = result_obj,
  assays   = assay_use,
  features = genes_use,
  group.by = "final_population",
  slot     = slot_use,
  verbose  = FALSE
)
avg_mat <- avg_list[[assay_use]]  # genes x populations

# percent detected per final_population (detected = expression > 0 in slot="data")
assay_use <- Seurat::DefaultAssay(result_obj)
slot_use  <- "data"

mat <- Seurat::GetAssayData(result_obj, assay = assay_use, slot = slot_use)[genes_use, , drop = FALSE]
grp <- result_obj$final_population

# split cells by population (drop NA explicitly)
idx_by_grp <- split(seq_along(grp), grp)
idx_by_grp <- idx_by_grp[!is.na(names(idx_by_grp))]

# percent detected per population (expression > 0 in slot="data")
pct_mat <- do.call(
  cbind,
  lapply(idx_by_grp, function(ii) Matrix::rowMeans(mat[, ii, drop = FALSE] > 0))
)

colnames(pct_mat) <- names(idx_by_grp)
rownames(pct_mat) <- rownames(mat)


colnames(avg_mat) <- sub("^g(?=\\p{Greek})", "", colnames(avg_mat), perl = TRUE)
grp_raw <- as.character(result_obj$final_population)
grp     <- sub("^g(?=\\p{Greek})", "", grp_raw, perl = TRUE)   # remove leading "g" only if next char is Greek
idx_by_grp <- split(seq_along(grp), grp)
idx_by_grp <- idx_by_grp[!is.na(names(idx_by_grp))]

# avg_mat and pct_mat built from the same grp names
mat <- Seurat::GetAssayData(result_obj, slot = "data")[genes_use, , drop = FALSE]

avg_mat <- do.call(cbind, lapply(idx_by_grp, function(ii) Matrix::rowMeans(mat[, ii, drop = FALSE])))
pct_mat <- do.call(cbind, lapply(idx_by_grp, function(ii) Matrix::rowMeans(mat[, ii, drop = FALSE] > 0)))

colnames(avg_mat) <- names(idx_by_grp); rownames(avg_mat) <- rownames(mat)
colnames(pct_mat) <- names(idx_by_grp); rownames(pct_mat) <- rownames(mat)

stopifnot(identical(colnames(avg_mat), colnames(pct_mat)))


# per-gene top population (top1, top2, delta)
top_tbl <- lapply(genes_use, function(g) {
  v <- avg_mat[g, ]
  ord <- order(v, decreasing = TRUE, na.last = TRUE)
  
  top1 <- names(v)[ord[1]]
  top2 <- if (length(ord) >= 2) names(v)[ord[2]] else NA_character_
  
  data.frame(
    gene      = g,
    top1_pop  = top1,
    top1_avg  = unname(v[ord[1]]),
    top1_pct  = unname(pct_mat[g, top1]),
    top2_pop  = top2,
    top2_avg  = if (!is.na(top2)) unname(v[ord[2]]) else NA_real_,
    top2_pct  = if (!is.na(top2)) unname(pct_mat[g, top2]) else NA_real_,
    delta12   = if (!is.na(top2)) unname(v[ord[1]] - v[ord[2]]) else NA_real_,
    stringsAsFactors = FALSE
  )
})

top_tbl <- dplyr::bind_rows(top_tbl) |>
  dplyr::arrange(dplyr::desc(delta12), dplyr::desc(top1_avg))

top_tbl
table(top_tbl$top1_pop)


Seurat::Idents(result_obj) <- "final_population"

feat_show <- top_tbl$gene[1:30]

?DotPlot
p <- Seurat::DotPlot(
  object   = result_obj,
  features = feat_show,
  group.by = "final_population",
  assay    = assay_use
) + ggplot2::coord_flip()+
  theme(axis.text.x = element_text(angle=70))

p


Seurat::Idents(result_obj) <- "final_population"
pops <- levels(Seurat::Idents(result_obj))


?FindMarkers
mk_list <- lapply(pops, function(pop) {
  m <- Seurat::FindMarkers(
    object   = result_obj,
    ident.1  = pop,
    ident.2  = NULL,          # vs all other populations
    features = genes_use,
    assay    = assay_use,
    slot     = slot_use,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    min.pct = 0.1,
  )
  m$gene <- rownames(m)
  m$pop  <- pop
  m
})

mk_tbl <- dplyr::bind_rows(mk_list)
mk_tbl <- mk_tbl[mk_tbl$p_val_adj<0.05,]


unique(mk_tbl$pop)


mk_tbl[mk_tbl$gene %in% genes_use,]
mk_tbl[mk_tbl$pct.1>0.5 & mk_tbl$pct.2<0.1,]








library(openxlsx)
cell_markers <- read.xlsx("G:/PhD_final/tables/cell_markers_curated_new_new_new_new.xlsx")
head(cell_markers)
result_obj
unique(result_obj$cluster_key_final)
unique(cell_markers[,c("Cell_population_general","Cell_population_detailed")])
# ---- Inputs ----
library(Seurat)
library(openxlsx)
library(dplyr)

obj <- result_obj
markers_xlsx <- "G:/PhD_final/tables/cell_markers_curated_new_new_new_new.xlsx"
cell_markers <- openxlsx::read.xlsx(markers_xlsx)

# Keep exact IDs, including ".1"; only remove whitespace/tabs
marker_raw <- unique(cell_markers$Markers_positive_SMESG)
marker_raw <- marker_raw[!is.na(marker_raw) & marker_raw != ""]
marker_raw <- trimws(gsub("\t", "", marker_raw))

assay_use <- "SCT"     # or "RNA"
layer_use <- "data"    # Seurat v5 layer

# Filter to markers present in the object (exact match, no trimming of ".1")
feat_all <- rownames(Seurat::GetAssayData(obj, assay = assay_use, layer = layer_use))
marker_present <- marker_raw[marker_raw %in% feat_all]
marker_missing <- setdiff(marker_raw, marker_present)

if (length(marker_missing) > 0) {
  message("Markers missing in ", assay_use, "/", layer_use, ": ", length(marker_missing))
}

# n_cells per cluster
n_cells_tbl <- as.data.frame(table(obj$cluster_key_final), stringsAsFactors = FALSE)
colnames(n_cells_tbl) <- c("cluster_key_final", "n_cells")

# mean_expr + pct_expr per cluster via DotPlot
dp <- Seurat::DotPlot(
  object   = obj,
  features = marker_present,
  group.by = "cluster_key_final",
  assay    = assay_use
)$data

# IMPORTANT: repair duplicate column names (e.g., duplicated UniProt.gene.symbol)
names(cell_markers) <- make.unique(names(cell_markers))

# Keep only the join-relevant columns (avoids carrying messy extras)
cell_markers_join <- cell_markers[, c(
  "Markers_positive_SMESG",
  "Cell_population_general",
  "Cell_population_detailed",
  "Markers_positive_common.name"
)]

# Clean only the marker IDs (keep ".1" etc intact)
cell_markers_join$Markers_positive_SMESG <- trimws(gsub("\t", "", cell_markers_join$Markers_positive_SMESG))

# Drop exact duplicate rows if they exist
cell_markers_join <- unique(cell_markers_join)

# Now the pipeline will work
stats_tbl <- dp %>%
  dplyr::transmute(
    cluster_key_final = .data$id,
    gene              = .data$features.plot,
    mean_expr         = .data$avg.exp,
    pct_expr          = .data$pct.exp
  ) %>%
  dplyr::left_join(n_cells_tbl, by = "cluster_key_final") %>%
  dplyr::left_join(cell_markers_join, by = c("gene" = "Markers_positive_SMESG")) %>%
  dplyr::arrange(cluster_key_final, dplyr::desc(pct_expr), dplyr::desc(mean_expr))

openxlsx::write.xlsx(
  stats_tbl,
  file = "G:/PhD_final/tables/cluster_marker_metrics_by_cluster_key_final.xlsx",
  overwrite = TRUE
)
head(stats_tbl)



library(dplyr)

thr <- 20  


cluster_marker_summary <- stats_tbl %>%
  dplyr::mutate(
    Markers_positive_common.name = trimws(gsub("\\s+", " ", Markers_positive_common.name))
  ) %>%
  dplyr::filter(
    !is.na(Markers_positive_common.name),
    Markers_positive_common.name != "",
    pct_expr >= thr
  ) %>%
  dplyr::group_by(cluster_key_final) %>%
  # disambiguate duplicated common names within each cluster by adding the SMESG ID
  dplyr::mutate(
    name_n_in_cluster = ave(Markers_positive_common.name, Markers_positive_common.name, FUN = length),
    marker_label = ifelse(
      name_n_in_cluster > 1,
      paste0(Markers_positive_common.name, " [", gene, "]"),
      Markers_positive_common.name
    )
  ) %>%
  dplyr::summarise(
    n_cells = first(n_cells),
    n_markers_passing = n_distinct(Markers_positive_common.name),
    markers_common = paste(
      marker_label[order(-pct_expr, marker_label)],
      sprintf("(%.1f%%)", pct_expr[order(-pct_expr, marker_label)]),
      sep = " ",
      collapse = "; "
    ),
    markers_SMESG = paste(sort(unique(gene)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_markers_passing), cluster_key_final)


head(cluster_marker_summary, 20)

as.data.frame(cluster_marker_summary)
head(cluster_marker_summary)

#unique(obj@meta.data$cluster_key_final,obj@meta.data$final_population_fixed)
cluster_marker_summary


library(dplyr)

md <- obj@meta.data %>%
  dplyr::select(cluster_key_final, final_population_fixed) %>%
  dplyr::filter(!is.na(cluster_key_final), !is.na(final_population_fixed))

# Mode + purity per cluster_key_final
cluster_old_anno <- md %>%
  dplyr::count(cluster_key_final, final_population_fixed, name = "n_label") %>%
  dplyr::group_by(cluster_key_final) %>%
  dplyr::mutate(
    n_cluster = sum(n_label),
    purity = n_label / n_cluster
  ) %>%
  dplyr::slice_max(order_by = n_label, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    cluster_key_final,
    final_population_fixed_mode = final_population_fixed,
    final_population_fixed_purity = round(100 * purity, 1),
    n_cells_meta = n_cluster
  )

# Join to your summary
cluster_marker_summary2 <- cluster_marker_summary %>%
  dplyr::left_join(cluster_old_anno, by = "cluster_key_final") %>%
  dplyr::relocate(final_population_fixed_mode, final_population_fixed_purity, .after = cluster_key_final) %>%
  dplyr::arrange(desc(n_markers_passing), cluster_key_final)

head(cluster_marker_summary2, 20)
openxlsx::write.xlsx(
  cluster_marker_summary2,
  file = "G:/PhD_final/tables/cluster_marker_summary2.xlsx",
  overwrite = TRUE
)










################################################################################

curated_anno <- openxlsx::read.xlsx("G:/PhD_final/tables/cluster_marker_summary_verified.xlsx")
head(curated_anno)
result_obj
head(result_obj$cluster_key_final)
result_obj
# map curated_anno$new_anno onto cells by their cluster label
cluster_key <- "cluster_key_final"   # change if you used a different column

colnames(curated_anno)[2] <- "final_population"
colnames(curated_anno)[1] <- "cluster"
stopifnot(all(c("cluster","final_population") %in% names(curated_anno)))
stopifnot(cluster_key %in% colnames(result_obj@meta.data))

# one row per cluster (in case of accidental duplicates)
cur_map_df <- curated_anno[!duplicated(curated_anno$cluster),
                           c("cluster","final_population")]
map <- setNames(as.character(cur_map_df$final_population),
                as.character(cur_map_df$cluster))

# cluster labels per cell
cl <- as.character(result_obj[[cluster_key, drop = TRUE]])

# annotate each cell
anno <- unname(map[cl])
anno[is.na(anno)] <- "Unknown"   # optional fallback

# store in metadata
result_obj$final_population <- factor(anno, levels = unique(cur_map_df$final_population))

# (optional) make these the current identities
# Seurat::Idents(result_obj) <- result_obj$curated_anno

# quick sanity check
table(result_obj$final_population, useNA = "ifany")
result_obj$final_population <- str_trim(result_obj$final_population)
save(result_obj,file="G:/PhD_final/result_obj_new.RData")



#cluster_marker_summary2_df <- as.data.frame(cluster_marker_summary2)
#cluster_marker_summary2_df[c(1:20,)]
# thresholds <- c(20, 30, 40, 50, 60, 70, 80)
# 
# threshold_scan <- lapply(thresholds, function(thr) {
#   stats_tbl %>%
#     dplyr::filter(!is.na(Markers_positive_common.name), Markers_positive_common.name != "", pct_expr >= thr) %>%
#     dplyr::group_by(cluster_key_final) %>%
#     dplyr::summarise(n_markers_passing = n_distinct(Markers_positive_common.name), .groups = "drop") %>%
#     dplyr::mutate(threshold = thr)
# }) %>%
#   dplyr::bind_rows() %>%
#   dplyr::arrange(threshold, desc(n_markers_passing))
# 
# threshold_scan %>% group_by(threshold) %>% summarise(
#   clusters_with_any = sum(n_markers_passing > 0),
#   median_markers = median(n_markers_passing),
#   .groups = "drop"
# )







# 
# # Option A: pick the single top cluster (no ties)
# top1 <- cluster_summ |> dplyr::slice_max(mean_expr, n = 1, with_ties = FALSE)
# abra_clusters_top1 <- top1$cluster_key_final
# 
# # Option B: pick a small set of "high" clusters (top 5 by mean)
# abra_clusters_top5 <- (cluster_summ |> dplyr::slice_max(mean_expr, n = 5, with_ties = TRUE))$cluster_key_final
# 
# abra_clusters_top1
# abra_clusters_top5
# # Choose which set you want to use:
# abra_clusters <- abra_clusters_top5
# # abra_clusters <- abra_clusters_top3
# 
# result_obj$final_population_fixed_abra <- as.character(result_obj$final_population_fixed)
# 
# idx <- !is.na(result_obj$cluster_key_final) & (result_obj$cluster_key_final %in% abra_clusters)
# result_obj$final_population_fixed_abra[idx] <- "Abraçada cell"
# 
# table(before = result_obj$final_population_fixed, after = result_obj$final_population_fixed_abra)[1:10, 1:10]
# table(result_obj$final_population_fixed_abra, useNA = "ifany")
# result_obj$final_population_fixed <- result_obj$final_population_fixed_abra
#save(result_obj,file="G:/PhD_final/result_obj_new.RData") 
# Seurat::VlnPlot(
#   result_obj,
#   features = gene_marker,
#   group.by = "cluster_key_final",
#   pt.size = 0
# )
# Seurat::VlnPlot(
#   result_obj,
#   features = gene_marker,
#   group.by = "final_population_fixed_abra",
#   pt.size = 0
# )
FAMILY_COLORS <- c(
  "σ-neoblast"                 = "#E3E3E3",
  "ν-neoblast"                 = "#747474",
  "ζ-neoblast"                 = "#848484",
  "γ-neoblast"                 = "#949494",
  "Eye neoblast"               = "#B9A5DD",
  "Muscle neoblast"            = "#C3C3C3",
  "Parenchymal neoblast"       = "#B3B3B3",
  "Protonephridial neoblast"   = "#A3A3A3",
  "Epidermal progenitor"       = "#ADCFE0",
  "Eye progenitor"             = "#00D3D3",
  "Pharyngeal progenitor"      = "#CD5D5E",
  "Pharyngeal neoblast"        = "#7F7F7F",
  "Parenchymal progenitor"     = "#F9B45C",
  "Protonephridial progenitor" = "#876173",
  "Epidermis"                  = "#0040FF",
  "Eye"                        = "#B9A5DD",
  "Nervous system"             = "#A084CA",
  "Pharynx"                    = "#D47072",
  "Muscle"                     = "#B71C1C",
  "Protonephridia"             = "#874A68",
  "Intestine"                  = "#3F783E",
  "Parenchyma"                 = "#F7A74C"
)
names(FAMILY_COLORS)[order(names(FAMILY_COLORS))]
# DETAILED_COLS <- c(
#   "Early epidermal progenitor"            = "#ABCEDF",
#   "Epidermal secretory gland"             = "#99BFD8",
#   "Epidermal secretory gland progenitor"  = "#88B0D1",
#   "Injury-induced epidermis"              = "#77A2CA",
#   "Late epidermal progenitor"             = "#6693C4",
#   "Mature epidermis (broad)"              = "#5685BD",
#   "Mature epidermis (DV-boundary)"        = "#4576B6",
#   "Mature epidermis (multiciliated)"      = "#3468B0",
#   "Eye photoreceptor"                     = "#00FFFF",
#   "Eye progenitor"                        = "#00D3D3",
#   "Pigment cup cell"                      = "#00A8A8",
#   "Basal cell"                            = "#D3E4BA",
#   "Goblet cell"                           = "#8DB180",
#   "Phagocyte"                             = "#497F46",
#   "Abraçada cell"                         = "#003c00",
#   "Anterior pole cell"                    = "#FFEBEE",
#   "BWM (circular)"                        = "#F2C7CA",
#   "BWM (dorsal midline )"                 = "#E6A3A6",
#   "BWM (longitudinal)"                    = "#D97F82",
#   "ECM-producing muscle"                  = "#CC5C5D",
#   "Posterior pole/PCG muscle"             = "#C1393A",
#   "Eye lineage neoblast"                  = "#838B8B",
#   "Muscle neoblast"                       = "#D5D5D5",
#   "Parenchymal neoblast"                  = "#C7C7C7",
#   "Protonephridial neoblast"              = "#B8B8B8",
#   "γ-neoblast (intestinal-fated)"         = "#AAAAAA",
#   "ζ-neoblast (dorsal epidermal-fated)"   = "#9B9B9B",
#   "ζ-neoblast (epidermal-fated)"          = "#8D8D8D",
#   "Pharyngeal neoblast"                   = "#7F7F7F",
#   "ν-neoblast (neural-fated)"             = "#707070",
#   "σ-neoblast (broad-lineage)"            = "#E3E3E3",
#   "Brain branch neuron"                   = "#D8CCF3",
#   "Catecholaminergic neuron"              = "#CDBFEB",
#   "Cholinergic neuron"                    = "#C3B2E4",
#   "GABAergic neuron"                      = "#B9A5DC",
#   "Glia"                                  = "#AF98D5",
#   "Glutamatergic neuron"                  = "#A58ACE",
#   "Glycinergic neuron"                    = "#9B7EC6",
#   "Mechanosensory neuron"                 = "#9171BF",
#   "Neuropeptidergic neuron"               = "#8764B8",
#   "OTF⁺ neuron"                           = "#7D57B1",
#   "PKD⁺ sensory neuron (ciliated)"        = "#7349A9",
#   "Serotonergic neuron"                   = "#693CA2",
#   "Tyraminergic neuron"                   = "#5F309B",
#   "AQP⁺ parenchymal cell"                 = "#F9E29D",
#   "LDLRR-1⁺ parenchymal cell"             = "#F3D38E",
#   "PGRN⁺ parenchymal cell"                = "#EFC57F",
#   "PSAP⁺ parenchymal cell"                = "#EBB670",
#   "PTF⁺ head parenchymal progenitor"      = "#E6A860",
#   "SSPO⁺ parenchymal cell"                = "#E19A51",
#   "SSPO⁺ parenchymal progenitor"          = "#DD8B42",
#   "Pharyngeal epithelium"                 = "#FFFF00",
#   "Pharyngeal muscle"                     = "#EBEB00",
#   "Pharyngeal progenitor"                 = "#D9D900",
#   "Protonephridial distal tubule cell"    = "#876173",
#   "Protonephridial flame cell"            = "#874A68",
#   "Protonephridial tubule cell"           = "#87345F",
#   "Protonephridial tubule precursor"      = "#871E56"
# )

DETAILED_COLS <- c(
  # --- Epidermis (blue) ---
  "Early epidermal progenitor"             = "#ABCEDF",
  "Late epidermal progenitor"              = "#6693C4",
  "Epidermis (broad)"                      = "#5685BD",
  "Epidermis (multiciliated)"              = "#3468B0",
  "Epidermal progenitor (multiciliated)"   = "#4576B6",
  "Epidermal secretory gland progenitor"   = "#99BFD8",
  
  # --- Eye (cyan) ---
  "Eye progenitor"                         = "#00D3D3",
  
  # --- Intestine / immune-like (green) ---
  "Basal cell"                             = "#D3E4BA",
  "Goblet cell"                            = "#8DB180",
  "Phagocyte (broad)"                     = "#497F46",
  
  
  # --- Pigment (brown) ---
  "Body pigment progenitor"                = "#C08A5A",
  "Body pigment cell"                      = "#8E5A3C",
  
  # --- Muscle (red) ---
  "Muscle progenitor"                      = "#F2C7CA",
  "BWM (dorsal midline)"                   = "#E6A3A6",
  "ECM-producing muscle"                   = "#CC5C5D",
  "Posterior pole/PCG muscle"              = "#C1393A",
  
  # --- Neoblasts (grey) ---
  "σ-neoblast (broad-lineage)"             = "#E3E3E3",
  "Muscle neoblast"                        = "#D5D5D5",
  "Protonephridial neoblast"               = "#B8B8B8",
  "Pharyngeal neoblast"                    = "#7F7F7F",
  "γ-neoblast (intestinal-fated)"          = "#AAAAAA",
  "ζ-neoblast (epidermal-fated)"           = "#8D8D8D",
  "ν-neoblast (neural-fated)"              = "#707070",
  "GLIRP-1⁺ parenchymal neoblast"          = "#B3B3B3",
  "PGRN⁺ parenchymal neoblast"             = "#BEBEBE",
  "FER3L-2⁺ parenchymal neoblast"          = "#C7C7C7",
  
  # --- Neural lineage (violet) ---
  "Neural progenitor (broad)"              = "#E9E4FA",
  "Glutamatergic neural progenitor"        = "#DED5F6",
  "Neuropeptidergic neural progenitor"     = "#D3C7F2",
  "Mechanosensory neural progenitor"       = "#C9BAEE",
  "PKD⁺ sensory neural progenitor"         = "#BFADE8",
  "Glia"                                   = "#AF98D5",
  "Brain branch neuron"                    = "#D8CCF3",
  "Catecholaminergic neuron"               = "#CDBFEB",
  "Cholinergic neuron"                     = "#C3B2E4",
  "Glutamatergic neuron"                   = "#A58ACE",
  "Mechanosensory neuron"                  = "#9171BF",
  "Neuropeptidergic neuron"                = "#8764B8",
  "PKD⁺ sensory neuron"                    = "#7349A9",
  "Serotonergic neuron"                    = "#693CA2",
  
  # --- Parenchyma (sand) ---
  "AQP⁺ parenchymal cell"                  = "#F9E29D",
  "LDLRR-1⁺ parenchymal cell"              = "#F3D38E",
  "GLIRP-1⁺ parenchymal progenitor"        = "#EDC281",
  "PSAP⁺ parenchymal progenitor"           = "#F1D6A0",
  "PSAP⁺ parenchymal cell"                 = "#EBB670",
  "PGRN⁺ parenchymal cell"                 = "#EFC57F",
  "FER3L-2⁺ parenchymal progenitor"        = "#E8B567",
  "NKX2⁺ parenchymal progenitor"           = "#E39F55",
  "PTF⁺ head parenchymal progenitor"       = "#E6A860",
  "SSPO⁺ parenchymal progenitor"           = "#E19A51",
  "SSPO⁺ parenchymal cell"                 = "#DD8B42",
  "Abraçada cell"                          = "#FAE8B4",
  
  # --- Pharynx (yellow) ---
  "Pharyngeal epithelium"                  = "#FFFF00",
  "Pharyngeal progenitor"                  = "#D9D900",
  "Pharyngeal phagocytic-type cell"        = "#C9C900",
  
  # --- Protonephridia (wine) ---
  "Protonephridial flame cell"             = "#874A68",
  "Protonephridial tubule cell"            = "#87345F"
)

# Optional sanity check:
stopifnot(setequal(unique(result_obj$final_population), names(DETAILED_COLS)))

names(DETAILED_COLS)[!(names(DETAILED_COLS) %in% unique(result_obj$final_population))]
unique(result_obj$final_population)[!(unique(result_obj$final_population) %in% names(DETAILED_COLS))]
#?DimPlot_scCustom
DimPlot_scCustom(seurat_object = result_obj,
                 colors_use = DETAILED_COLS,
                 reduction = "umap.d33.nn100.md0.3",
                 group.by = "final_population",
                 label = TRUE,
                 label.box = TRUE,
                 repel = TRUE,
                 pt.size = 1,
                 figure_plot = TRUE) + NoLegend()
unique(result_obj$final_population)
unique(result_obj$condition_pretty)



# =========================
# 100% stacked "condition-within-celltype" plot + hit outlining (no replicates)
# =========================

stopifnot(inherits(result_obj, "Seurat"))

md0 <- result_obj@meta.data
stopifnot(all(c("final_population", "condition_pretty") %in% colnames(md0)))

md <- tibble::as_tibble(md0, rownames = "cell")

has_cc <- "condition_correct" %in% colnames(md)

md <- md |>
  dplyr::mutate(
    final_population = as.character(.data$final_population),
    condition_pretty = as.character(.data$condition_pretty),
    
    condition_base = stringr::str_replace(.data$condition_pretty, "\\n.*$", ""),
    timepoint_h    = as.integer(stringr::str_remove(
      stringr::str_extract(.data$condition_pretty, "\\d+hpa"),
      "hpa"
    )),
    
    group = dplyr::case_when(
      .data$condition_base == "WT"       ~ "WT",
      .data$condition_base == "GFP Mock" ~ "GFP",
      .data$condition_base == "ELAC2 KD" ~ "ELAC",
      TRUE ~ NA_character_
    ),
    
    # build condition_code depending on whether condition_correct exists
    condition_code = if (has_cc) as.character(.data$condition_correct) else paste0(.data$group, .data$timepoint_h),
    
    #final_population_plot = stringr::str_replace_all(.data$final_population, "[^A-Za-z0-9]+", "_")
    final_population_plot = as.character(.data$final_population)
    
  ) |>
  dplyr::filter(
    !is.na(.data$group),
    !is.na(.data$timepoint_h),
    !is.na(.data$final_population_plot)
  )


# ---- settings for "descriptive" hit calling (no biological replicates) ----
alpha <- 0.05
min_abs_delta <- 0.02   # absolute prop change threshold (2 percentage points)
min_fold      <- 1.25   # fold threshold (set to 1 to disable)

two_prop_test <- function(x1, n1, x2, n2) {
  tab <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE)
  expected <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$expected)
  use_fisher <- any(expected < 5)
  
  p <- if (use_fisher) {
    stats::fisher.test(tab)$p.value
  } else {
    stats::prop.test(x = c(x1, x2), n = c(n1, n2), correct = FALSE)$p.value
  }
  
  p1 <- x1 / n1
  p2 <- x2 / n2
  eps <- 1e-9
  fc  <- (p1 + eps) / (p2 + eps)
  
  tibble::tibble(
    p_value = unname(p),
    prop_1 = p1, prop_2 = p2,
    delta_prop = p1 - p2,
    fold = fc
  )
}
COND_COLS <- c(
  WT0="#619CFF", GFP0="#00BA38", ELAC0="#F8766D",
  WT16="#619CFF", GFP16="#00BA38", ELAC16="#F8766D",
  WT24="#619CFF", GFP24="#00BA38", ELAC24="#F8766D",
  WT72="#619CFF", GFP72="#00BA38", ELAC72="#F8766D"
)


call_hits_one_timepoint <- function(md_tp) {
  cnt <- md_tp |>
    dplyr::count(.data$final_population_plot, .data$group, name = "x") |>
    tidyr::pivot_wider(names_from = .data$group, values_from = .data$x, values_fill = 0)
  
  # Ensure all three columns exist (pivot_wider drops missing groups)
  for (nm in c("WT", "GFP", "ELAC")) {
    if (!nm %in% colnames(cnt)) cnt[[nm]] <- 0L
  }
  
  n_WT  <- sum(md_tp$group == "WT")
  n_GFP <- sum(md_tp$group == "GFP")
  n_EL  <- sum(md_tp$group == "ELAC")
  stopifnot(n_WT > 0, n_GFP > 0, n_EL > 0)
  
  res <- cnt |>
    dplyr::rowwise() |>
    dplyr::mutate(
      tmp_el_wt  = list(two_prop_test(ELAC, n_EL, WT,  n_WT)),
      tmp_gfp_wt = list(two_prop_test(GFP,  n_GFP, WT,  n_WT)),
      tmp_el_gfp = list(two_prop_test(ELAC, n_EL, GFP, n_GFP))
    ) |>
    dplyr::ungroup() |>
    tidyr::unnest_wider(.data$tmp_el_wt,  names_sep = "_") |>
    dplyr::rename_with(~ sub("^tmp_el_wt_",  "ELACvsWT_",  .x), dplyr::starts_with("tmp_el_wt_")) |>
    tidyr::unnest_wider(.data$tmp_gfp_wt, names_sep = "_") |>
    dplyr::rename_with(~ sub("^tmp_gfp_wt_", "GFPvsWT_",   .x), dplyr::starts_with("tmp_gfp_wt_")) |>
    tidyr::unnest_wider(.data$tmp_el_gfp, names_sep = "_") |>
    dplyr::rename_with(~ sub("^tmp_el_gfp_", "ELACvsGFP_", .x), dplyr::starts_with("tmp_el_gfp_"))
  
  res <- res |>
    dplyr::mutate(
      fdr_ELACvsWT  = stats::p.adjust(.data$ELACvsWT_p_value,  method = "BH"),
      fdr_GFPvsWT   = stats::p.adjust(.data$GFPvsWT_p_value,   method = "BH"),
      fdr_ELACvsGFP = stats::p.adjust(.data$ELACvsGFP_p_value, method = "BH"),
      
      injection_affected = .data$fdr_GFPvsWT < alpha,
      
      same_dir_ELAC_WT_vs_ELAC_GFP =
        (sign(.data$ELACvsWT_delta_prop) == sign(.data$ELACvsGFP_delta_prop)) &
        (sign(.data$ELACvsWT_delta_prop) != 0),
      
      # NEW injection-correction logic
      elac2_specific_raw =
        (.data$fdr_ELACvsWT < alpha) &
        (
          (!.data$injection_affected) |
            (
              .data$injection_affected &
                (.data$fdr_ELACvsGFP < alpha) &
                .data$same_dir_ELAC_WT_vs_ELAC_GFP
            )
        ),
      
      pass_effect =
        (abs(.data$ELACvsWT_delta_prop) >= min_abs_delta) &
        (.data$ELACvsWT_fold >= min_fold | .data$ELACvsWT_fold <= (1 / min_fold)),
      
      elac2_specific = .data$elac2_specific_raw & .data$pass_effect
    )
  
  
  res
}


# ---- build plotting table: within each cell type, fractions sum to 1 (100% stacked) ----
df_stack <- md |>
  dplyr::count(.data$timepoint_h, .data$final_population_plot, .data$condition_code, name = "n") |>
  dplyr::group_by(.data$timepoint_h, .data$final_population_plot) |>
  dplyr::mutate(frac = .data$n / sum(.data$n)) |>
  dplyr::ungroup()

# keep an ordering similar to your PDF (alphabetical)
pop_levels <- sort(unique(df_stack$final_population_plot))
df_stack$final_population_plot <- factor(df_stack$final_population_plot, levels = pop_levels)

# ---- call hits per timepoint ----
timepoints <- sort(unique(df_stack$timepoint_h))
colnames(md)
hit_map <- lapply(timepoints, function(tp) {
  md_tp <- md |> dplyr::filter(.data$timepoint_h == tp)
  ht <- call_hits_one_timepoint(md_tp) |>
    dplyr::filter(.data$elac2_specific) |>
    dplyr::transmute(timepoint_h = tp, final_population_plot = .data$final_population_plot)
  ht
}) |>
  dplyr::bind_rows()

# rects need numeric x positions
rect_df <- hit_map |>
  dplyr::mutate(
    x = as.numeric(factor(.data$final_population_plot, levels = pop_levels)),
    xmin = .data$x - 0.5,
    xmax = .data$x + 0.5,
    ymin = 0,
    ymax = 1
  )
# mark which (timepoint, celltype) should be highlighted
df_stack2 <- df_stack |>
  dplyr::mutate(
    highlight = dplyr::if_else(
      paste(timepoint_h, final_population_plot) %in%
        paste(rect_df$timepoint_h, rect_df$final_population_plot),
      "hit", "other"
    )
  )

# ---- final plot (faceted like your example) ----
p <- ggplot2::ggplot(df_stack2, ggplot2::aes(
  x = .data$final_population_plot,
  y = .data$frac,
  fill = .data$condition_code,
  alpha = .data$highlight
)) +
  ggplot2::geom_col(width = 0.9) +
  ggplot2::facet_wrap(~ timepoint_h, ncol = 1, labeller = ggplot2::labeller(
    timepoint_h = function(x) paste0(x, " hpa")
  )) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  ggplot2::scale_alpha_manual(values = c(hit = 1, other = 0.25), guide = "none") +
  ggplot2::labs(
    x = NULL,
    y = "Relative cell type fraction",
    fill = "Condition"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0, size = 9),
    strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
    strip.text = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::scale_fill_manual(values = COND_COLS)

as.data.frame(df_stack2[df_stack2$highlight!="other",])
# outline highlighted cell types (injection-corrected ELAC2 vs WT hits)
if (nrow(rect_df) > 0) {
  p <- p +
    ggplot2::geom_rect(
      data = rect_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.4
    )
}

# apply your palette if available
if (exists("DETAILED_COLS", inherits = TRUE)) {
  p <- p + ggplot2::scale_fill_manual(values = DETAILED_COLS)
}

p <- p + ggplot2::scale_fill_manual(values = COND_COLS)

print(p)

# optionally:
# ggplot2::ggsave("cell_barplot_sig_like_example.pdf", p, width = 14, height = 9)

# save all plots to PDF
# grDevices::pdf("celltype_proportions_no_reps_ELAC2_injection_corrected.pdf", width = 11, height = 14)
# for (nm in names(plots)) if (!is.null(plots[[nm]])) print(plots[[nm]])
# grDevices::dev.off()


tp  <- 72
pop <- "Abraçada cell"

md_tp <- md |>
  dplyr::filter(.data$timepoint_h == tp)

# totals per group (within tp)
tot <- md_tp |>
  dplyr::count(.data$group, name = "n_total")

# counts for the chosen population (within tp)
cnt <- md_tp |>
  dplyr::filter(.data$final_population == pop) |>
  dplyr::count(.data$group, name = "x")

tab <- tot |>
  dplyr::left_join(cnt, by = "group") |>
  dplyr::mutate(
    x = dplyr::coalesce(.data$x, 0L),
    prop_within_condition = .data$x / .data$n_total
  ) |>
  dplyr::mutate(
    frac_within_celltype = .data$x / sum(.data$x)
  )

tab
x_WT  <- tab$x[tab$group == "WT"];   n_WT  <- tab$n_total[tab$group == "WT"]
x_GFP <- tab$x[tab$group == "GFP"];  n_GFP <- tab$n_total[tab$group == "GFP"]
x_EL  <- tab$x[tab$group == "ELAC"]; n_EL  <- tab$n_total[tab$group == "ELAC"]

t_el_wt  <- two_prop_test(x_EL,  n_EL,  x_WT,  n_WT)
t_gfp_wt <- two_prop_test(x_GFP, n_GFP, x_WT,  n_WT)
t_el_gfp <- two_prop_test(x_EL,  n_EL,  x_GFP, n_GFP)

list(
  ELAC_vs_WT  = t_el_wt,
  GFP_vs_WT   = t_gfp_wt,
  ELAC_vs_GFP = t_el_gfp
)

call_hits_one_timepoint(md_tp) |>
  dplyr::filter(.data$final_population_plot == pop) |>
  dplyr::select(
    final_population_plot, WT, GFP, ELAC,
    ELACvsWT_p_value, ELACvsWT_delta_prop, ELACvsWT_fold, fdr_ELACvsWT,
    GFPvsWT_p_value,  GFPvsWT_delta_prop,  GFPvsWT_fold,  fdr_GFPvsWT,
    elac2_specific_raw, pass_effect, elac2_specific
  )









#################################################################################



meta(result_obj)
result_obj
meta <- result_obj@meta.data
colnames(meta)
unique(result_obj$condition_pretty)
unique(result_obj$condition_correct)
library(dplyr)
library(ggplot2)
library(forcats)
library(cowplot)

# ---- palette (your input) ----
COND_COLS <- c(
  "WT\n0hpa"       = "#619CFF",
  "GFP Mock\n0hpa" = "#00BA38",
  "ELAC2 KD\n0hpa" = "#F8766D",
  "WT\n16hpa"       = "#619CFF",
  "GFP Mock\n16hpa" = "#00BA38",
  "ELAC2 KD\n16hpa" = "#F8766D",
  "WT\n24hpa"       = "#619CFF",
  "GFP Mock\n24hpa" = "#00BA38",
  "ELAC2 KD\n24hpa" = "#F8766D",
  "WT\n72hpa"       = "#619CFF",
  "GFP Mock\n72hpa" = "#00BA38",
  "ELAC2 KD\n72hpa" = "#F8766D"
)

# Replace your cond_plain() with this version
cond_plain <- function(x) {
  x <- gsub("\n", " ", x, fixed = TRUE)                 # WT\n0hpa -> WT 0hpa
  x <- gsub("\\bMock\\b", "mock", x)                    # Mock -> mock (word-boundary safe)
  x <- gsub("([0-9]+)hpa\\b", "\\1 hpa", x, perl = TRUE) # 72hpa -> 72 hpa
  x <- sub("^ELAC2 KD", "Smed ELAC2 KD", x)             # add Smed
  x
}

# If you are using plotmath italics, keep this (no change needed except it now uses updated cond_plain)
cond_axis_expr <- function(x) {
  x <- cond_plain(x)
  if (grepl("^Smed ELAC2", x)) {
    suffix <- sub("^Smed ELAC2", "", x)  # e.g. " KD 72 hpa"
    paste0('paste(italic("Smed ELAC2"), "', suffix, '")')
  } else {
    paste0('"', x, '"')
  }
}

# Rebuild palette names to match the transformed labels
COND_COLS2 <- COND_COLS
names(COND_COLS2) <- cond_plain(names(COND_COLS2))
# ---- rebuild df with matching factor levels ----
df <- result_obj@meta.data %>%
  mutate(
    condition_pretty = factor(
      cond_plain(as.character(condition_pretty)),
      levels = names(COND_COLS2)
    ),
    nCount_RNA = as.numeric(nCount_RNA),
    nFeature_RNA = as.numeric(nFeature_RNA),
    percent.mt = as.numeric(percent.mt),
    log10GenesPerUMI = as.numeric(log10GenesPerUMI),
    log10_nCount = log10(nCount_RNA + 1),
    log10_nFeature = log10(nFeature_RNA + 1),
    umi_per_gene = nCount_RNA / pmax(nFeature_RNA, 1),
    log10_umi_per_gene = log10(umi_per_gene + 1e-6)
  )
# Rebuild df factor levels to match COND_COLS2 exactly
df <- df %>%
  mutate(
    condition_pretty = factor(
      cond_plain(as.character(condition_pretty)),
      levels = names(COND_COLS2)
    )
  )
stopifnot(identical(levels(df$condition_pretty), names(COND_COLS2)))

# ---- theme (no ggtext needed) ----
theme_qc <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(face = "bold")
  )

# ---- helper for pretty violins ----
vln_pretty <- function(data, y, ylab) {
  ggplot(data, aes(x = condition_pretty, y = .data[[y]], fill = condition_pretty)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.95) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, alpha = 0.25, color = "grey35") +
    stat_summary(fun = median, geom = "point", size = 1.0, color = "black") +
    coord_flip() +
    scale_fill_manual(values = COND_COLS2, drop = FALSE) +
    scale_x_discrete(labels = function(x) parse(text = vapply(x, cond_axis_expr, character(1)))) +
    labs(y = ylab) +
    theme_qc
}

# ---- cells retained ----
p_cells <- df %>%
  count(condition_pretty, name = "n_cells") %>%
  ggplot(aes(x = condition_pretty, y = n_cells, fill = condition_pretty)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = n_cells), hjust = -0.05, size = 3.0) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = COND_COLS2, drop = FALSE) +
  scale_x_discrete(labels = function(x) parse(text = vapply(x, cond_axis_expr, character(1)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(y = "Cells retained") +
  theme_qc

# ---- other QC panels ----
p_umi   <- vln_pretty(df, "log10_nCount",       "Library size: log10(UMI + 1)")
p_genes <- vln_pretty(df, "log10_nFeature",     "Detected genes: log10(genes + 1)")
p_mt    <- vln_pretty(df, "percent.mt",         "Mitochondrial reads (%)")
p_gpu   <- vln_pretty(df, "log10GenesPerUMI",   "Complexity: log10(genes per UMI)")

# Add a bit more top margin to every plot (so labels don't hit y tick labels)
add_top_margin <- function(p, top = 10) {
  p + theme(plot.margin = margin(t = top, r = 6, b = 6, l = 6))
}

p_cells2 <- add_top_margin(p_cells, 12)
p_umi2   <- add_top_margin(p_umi,   12)
p_genes2 <- add_top_margin(p_genes, 12)
p_mt2    <- add_top_margin(p_mt,    12)

p_all <- cowplot::plot_grid(
  p_cells2, p_umi2,
  p_genes2, p_mt2,
  ncol = 2,
  align = "hv",
  labels = c("A","B","C","D"),
  label_size = 12,
  label_x = 0.02,   # a little right
  label_y = 1.01    # inside the panel
)

p_all


# ggsave("QC_by_condition_pretty_6panel.pdf", p_all, width = 12, height = 10)








# ggsave("QC_by_condition_pretty_6panel.pdf", p_all, width = 12, height = 10)

# ggsave("QC_by_condition_pretty.pdf", p_all, width = 9, height = 14)


#################################################################################
result_obj
table(result_obj$final_population)

options(ggrepel.max.overlaps = Inf)
p <- DimPlot_scCustom(
  seurat_object = result_obj,
  colors_use = DETAILED_COLS,
  reduction = "umap.d33.nn100.md0.3",
  group.by = "final_population",
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 0.5,
  figure_plot = TRUE
)

p & ggplot2::theme(legend.position = "none")
red <- "umap.d33.nn100.md0.3"

emb <- as.data.frame(Seurat::Embeddings(result_obj, reduction = red))
colnames(emb)
# typically two columns, e.g. "UMAP_1" "UMAP_2" (or similar)

x <- emb[[1]]
y <- emb[[2]]

bad_idx <- which(!is.finite(x) | !is.finite(y) | is.na(x) | is.na(y))

length(bad_idx)
bad_cells <- rownames(emb)[bad_idx]
bad_cells

emb[bad_idx, , drop = FALSE]
result_obj@meta.data[bad_cells, c("final_population"), drop = FALSE]

# Seurat::FetchData(result_obj, vars = c(colnames(emb)[1], colnames(emb)[2], "final_population")) |>
#   dplyr::filter(!is.finite(.data[[1]]) | !is.finite(.data[[2]]) | is.na(.data[[1]]) | is.na(.data[[2]]) | is.na(final_population))
df <- Seurat::FetchData(
  result_obj,
  vars = c(colnames(emb)[1], colnames(emb)[2], "final_population")
)

c1 <- colnames(emb)[1]
c2 <- colnames(emb)[2]

bad <- df |>
  dplyr::filter(
    !is.finite(.data[[c1]]) | !is.finite(.data[[c2]]) |
      is.na(.data[[c1]]) | is.na(.data[[c2]]) |
      is.na(final_population)
  )

nrow(bad)
rownames(bad)
bad



?FeaturePlot_scCustom
FeaturePlot_scCustom(seurat_object = result_obj,
                     features = "SMESG000049716.1",
                 reduction = "umap.d33.nn100.md0.3",
                 )



(our_annotation=rtracklayer::import('D:/Elac2/final_anno_with_mtDNA_and_ncRNA_new.gtf'))
our_annotation_tRNA <- our_annotation[our_annotation$annotation=="tRNA",]
unique(our_annotation_tRNA$detailed_annotation)
our_annotation_tRNA_Gly <- our_annotation_tRNA[grepl("Gly",our_annotation_tRNA$detailed_annotation) |
                                                 grepl("GLY",our_annotation_tRNA$detailed_annotation), ]
unique(our_annotation_tRNA_Gly$detailed_annotation)
our_annotation_tRNA_Gly[our_annotation_tRNA_Gly$detailed_annotation=="dd_Smed_g4_Pseudo_GLY-GCC_tRNA_11",]

