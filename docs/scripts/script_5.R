
# Marker statistics per cluster 
#   1) Computes average expression and % detected for a given gene set per cluster.
#   2) Builds a long, join-ready marker table from Seurat::DotPlot() data.
#   3) Builds a per-cluster marker summary (top markers passing a % threshold).
#   4) Optionally maps a curated per-cluster annotation file back to the Seurat object.

# ----------------------------
# 0) Minimal dependencies
# ----------------------------
stopifnot(requireNamespace("Seurat", quietly = TRUE))
stopifnot(requireNamespace("Matrix", quietly = TRUE))
stopifnot(requireNamespace("dplyr", quietly = TRUE))
stopifnot(requireNamespace("openxlsx", quietly = TRUE))

# ----------------------------
# 1) Small helpers
# ----------------------------

# Seurat v4/v5 compatibility: GetAssayData() uses slot (v4) or layer (v5).
.get_data <- function(obj, assay, layer_or_slot = "data") {
  # Seurat v5 uses `layer=`; v4 uses `slot=`.
  fmls <- names(formals(Seurat::GetAssayData))
  if ("layer" %in% fmls) {
    return(Seurat::GetAssayData(obj, assay = assay, layer = layer_or_slot))
  }
  Seurat::GetAssayData(obj, assay = assay, slot = layer_or_slot)
}

# Remove leading "g" only when it precedes a Greek letter (to keep g1/g2 etc intact).
.strip_g_before_greek <- function(x) {
  x <- as.character(x)
  sub("^g(?=\\p{Greek})", "", x, perl = TRUE)
}

# Safe column extraction from meta.data
.get_md <- function(obj, col) {
  stopifnot(col %in% colnames(obj@meta.data))
  obj@meta.data[[col]]
}

# ----------------------------
# 2) Core computations
# ----------------------------

# Compute average expression + percent detected per group for gene list.
#
compute_avg_and_pct <- function(obj,
                                genes,
                                group_col,
                                assay = Seurat::DefaultAssay(obj),
                                layer_or_slot = "data",
                                fix_group_names = TRUE) {
  stopifnot(is.character(genes), length(genes) > 0)
  stopifnot(group_col %in% colnames(obj@meta.data))

  genes_use <- intersect(genes, rownames(obj))
  if (!length(genes_use)) stop("None of the provided genes are present in the object.", call. = FALSE)

  grp_raw <- .get_md(obj, group_col)
  grp <- as.character(grp_raw)
  if (fix_group_names) grp <- .strip_g_before_greek(grp)

  # Drop NA groups 
  keep_cells <- !is.na(grp)
  grp <- grp[keep_cells]

  mat <- .get_data(obj, assay = assay, layer_or_slot = layer_or_slot)
  mat <- mat[genes_use, keep_cells, drop = FALSE]

  idx_by_grp <- split(seq_along(grp), grp)

  avg_mat <- do.call(cbind, lapply(idx_by_grp, function(ii) {
    Matrix::rowMeans(mat[, ii, drop = FALSE])
  }))

  pct_mat <- do.call(cbind, lapply(idx_by_grp, function(ii) {
    Matrix::rowMeans(mat[, ii, drop = FALSE] > 0)
  }))

  colnames(avg_mat) <- names(idx_by_grp)
  colnames(pct_mat) <- names(idx_by_grp)
  rownames(avg_mat) <- rownames(mat)
  rownames(pct_mat) <- rownames(mat)

  stopifnot(identical(colnames(avg_mat), colnames(pct_mat)))

  list(
    genes_use = genes_use,
    group_levels = colnames(avg_mat),
    avg_mat = avg_mat,
    pct_mat = pct_mat
  )
}

# For each gene: report top1/top2 groups and the delta between them.
rank_top_groups_per_gene <- function(avg_mat, pct_mat) {
  stopifnot(all(dim(avg_mat) == dim(pct_mat)))

  genes <- rownames(avg_mat)
  out <- lapply(genes, function(g) {
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

  dplyr::bind_rows(out) |>
    dplyr::arrange(dplyr::desc(.data$delta12), dplyr::desc(.data$top1_avg))
}

# Build a long marker table for a marker reference XLSX using DotPlot() output.
build_marker_stats_table <- function(obj,
                                     markers_xlsx,
                                     marker_gene_col = "Markers_positive_SMESG",
                                     group_col = "cluster_key_final",
                                     assay = "SCT",
                                     layer_or_slot = "data") {
  stopifnot(file.exists(markers_xlsx))
  stopifnot(group_col %in% colnames(obj@meta.data))

  # Read marker reference and keep only the relevant join columns.
  mk <- openxlsx::read.xlsx(markers_xlsx)
  names(mk) <- make.unique(names(mk))

  need_cols <- c(
    marker_gene_col,
    "Cell_population_general",
    "Cell_population_detailed",
    "Markers_positive_common.name"
  )
  miss <- setdiff(need_cols, names(mk))
  if (length(miss)) stop("Missing columns in marker XLSX: ", paste(miss, collapse = ", "), call. = FALSE)

  mk_join <- mk[, need_cols, drop = FALSE]
  mk_join[[marker_gene_col]] <- trimws(gsub("\\t", "", mk_join[[marker_gene_col]]))
  mk_join <- unique(mk_join)

  # Restrict to markers present in the object (exact match; keep ".1" etc intact).
  feat_all <- rownames(.get_data(obj, assay = assay, layer_or_slot = layer_or_slot))
  marker_raw <- unique(mk_join[[marker_gene_col]])
  marker_raw <- marker_raw[!is.na(marker_raw) & marker_raw != ""]
  marker_present <- marker_raw[marker_raw %in% feat_all]

  if (!length(marker_present)) stop("No marker genes from XLSX were found in the object (assay/layer mismatch?).", call. = FALSE)

  # Cluster size table
  n_cells_tbl <- as.data.frame(table(obj[[group_col, drop = TRUE]]), stringsAsFactors = FALSE)
  colnames(n_cells_tbl) <- c(group_col, "n_cells")

  # DotPlot gives per-group mean and % detected (in the plotting data). 
  dp <- Seurat::DotPlot(
    object   = obj,
    features = marker_present,
    group.by = group_col,
    assay    = assay
  )$data

  # Normalize column names across Seurat versions
  stats_tbl <- dp |>
    dplyr::transmute(
      !!group_col := .data$id,
      gene        = .data$features.plot,
      mean_expr   = .data$avg.exp,
      pct_expr    = .data$pct.exp
    ) |>
    dplyr::left_join(n_cells_tbl, by = group_col) |>
    dplyr::left_join(mk_join, by = setNames(marker_gene_col, "gene")) |>
    dplyr::arrange(.data[[group_col]], dplyr::desc(.data$pct_expr), dplyr::desc(.data$mean_expr))

  stats_tbl
}

# Summarize markers per cluster for quick manual review.
build_cluster_marker_summary <- function(stats_tbl,
                                         group_col = "cluster_key_final",
                                         thr = 20) {
  stopifnot(all(c(group_col, "gene", "pct_expr", "mean_expr", "n_cells") %in% colnames(stats_tbl)))

  stats_tbl |>
    dplyr::mutate(
      Markers_positive_common.name = trimws(gsub("\\s+", " ", .data$Markers_positive_common.name))
    ) |>
    dplyr::filter(
      !is.na(.data$Markers_positive_common.name),
      .data$Markers_positive_common.name != "",
      .data$pct_expr >= thr
    ) |>
    dplyr::group_by(.data[[group_col]]) |>
   
    dplyr::mutate(
      name_n_in_cluster = ave(.data$Markers_positive_common.name, .data$Markers_positive_common.name, FUN = length),
      marker_label = ifelse(
        .data$name_n_in_cluster > 1,
        paste0(.data$Markers_positive_common.name, " [", .data$gene, "]"),
        .data$Markers_positive_common.name
      )
    ) |>
    dplyr::summarise(
      n_cells = dplyr::first(.data$n_cells),
      n_markers_passing = dplyr::n_distinct(.data$Markers_positive_common.name),
      markers_common = paste(
        .data$marker_label[order(-.data$pct_expr, .data$marker_label)],
        sprintf("(%.1f%%)", .data$pct_expr[order(-.data$pct_expr, .data$marker_label)]),
        sep = " ",
        collapse = "; "
      ),
      markers_SMESG = paste(sort(unique(.data$gene)), collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_markers_passing), .data[[group_col]])
}

# Compute mode + purity of an existing annotation per cluster.
mode_and_purity_by_cluster <- function(obj, cluster_col, anno_col) {
  md <- obj@meta.data |>
    dplyr::select(dplyr::all_of(c(cluster_col, anno_col))) |>
    dplyr::filter(!is.na(.data[[cluster_col]]), !is.na(.data[[anno_col]]))

  md |>
    dplyr::count(.data[[cluster_col]], .data[[anno_col]], name = "n_label") |>
    dplyr::group_by(.data[[cluster_col]]) |>
    dplyr::mutate(
      n_cluster = sum(.data$n_label),
      purity = .data$n_label / .data$n_cluster
    ) |>
    dplyr::slice_max(order_by = .data$n_label, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      !!cluster_col := .data[[cluster_col]],
      anno_mode = .data[[anno_col]],
      anno_purity = round(100 * .data$purity, 1),
      n_cells_meta = .data$n_cluster
    )
}

# Map curated per-cluster annotations back onto cell
apply_curated_cluster_annotation <- function(obj,
                                             curated_xlsx,
                                             cluster_col = "cluster_key_final",
                                             out_col = "final_population",
                                             fallback = "Unknown") {
  stopifnot(file.exists(curated_xlsx))
  stopifnot(cluster_col %in% colnames(obj@meta.data))

  cur <- openxlsx::read.xlsx(curated_xlsx)
  if (ncol(cur) < 2) stop("Curated XLSX must have at least 2 columns: cluster and final_population.", call. = FALSE)

  # Standardize expected names: first col = cluster, second col = final_population.
  names(cur)[1:2] <- c("cluster", "final_population")
  cur <- cur[!duplicated(cur$cluster) & !is.na(cur$cluster), c("cluster", "final_population"), drop = FALSE]

  map <- setNames(as.character(cur$final_population), as.character(cur$cluster))
  cl <- as.character(obj[[cluster_col, drop = TRUE]])

  anno <- unname(map[cl])
  anno[is.na(anno) | anno == ""] <- fallback

  obj[[out_col]] <- factor(anno, levels = unique(cur$final_population))
  obj
}



# Parameters
cluster_col <- "cluster_key_final"      # cluster labels
assay_use   <- "SCT"                    # or Seurat::DefaultAssay(result_obj)
layer_use   <- "data"                   # "data" for log-normalized/SCT

# Marker XLSX workflow
markers_xlsx <- "G:/PhD_final/tables/cell_markers_curated_new_new_new_new.xlsx"

# Output paths
out_stats_xlsx   <- "G:/PhD_final/tables/cluster_marker_metrics_by_cluster_key_final.xlsx"
out_summary_xlsx <- "G:/PhD_final/tables/cluster_marker_summary2.xlsx"

# Optional curated mapping
curated_xlsx <- "G:/PhD_final/tables/cluster_marker_summary_verified.xlsx"
out_obj_rdata <- "G:/PhD_final/result_obj_new.RData"

# --- Build marker statistics and write to XLSX ---
stats_tbl <- build_marker_stats_table(
  obj = result_obj,
  markers_xlsx = markers_xlsx,
  marker_gene_col = "Markers_positive_SMESG",
  group_col = cluster_col,
  assay = assay_use,
  layer_or_slot = layer_use
)
openxlsx::write.xlsx(stats_tbl, file = out_stats_xlsx, overwrite = TRUE)

# --- Summarize markers per cluster (for manual review) ---
cluster_marker_summary <- build_cluster_marker_summary(stats_tbl, group_col = cluster_col, thr = 20)

anno_col_existing <- "final_population_fixed"
if (anno_col_existing %in% colnames(result_obj@meta.data)) {
  old_anno <- mode_and_purity_by_cluster(result_obj, cluster_col = cluster_col, anno_col = anno_col_existing)
  cluster_marker_summary <- cluster_marker_summary |>
    dplyr::left_join(old_anno, by = cluster_col) |>
    dplyr::relocate(.data$anno_mode, .data$anno_purity, .after = dplyr::all_of(cluster_col))
}

openxlsx::write.xlsx(cluster_marker_summary, file = out_summary_xlsx, overwrite = TRUE)

# --- Apply curated annotation back to Seurat object ---
if (file.exists(curated_xlsx)) {
  result_obj <- apply_curated_cluster_annotation(
    obj = result_obj,
    curated_xlsx = curated_xlsx,
    cluster_col = cluster_col,
    out_col = "final_population",
    fallback = "Unknown"
  )

  # Basic cleanup 
  result_obj$final_population <- factor(trimws(as.character(result_obj$final_population)))

  save(result_obj, file = out_obj_rdata)
}
