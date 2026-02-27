# Purpose
#   Evaluate integration quality on a chosen embedding (usually PCA):
#     - iLISI: local mixing by batch/condition (higher = better mixing)
#     - cLISI-derived purity: local purity by cell type (higher = better purity)
#     - kBET: average rejection rate of local batch-mixing tests (lower = better mixing)


# ----------------------------
# 0) Package checks
# ----------------------------
check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package: ", pkg, " (install it first).", call. = FALSE)
  }
}
check_pkg("Seurat")
check_pkg("lisi")   # immunogenomics/LISI (R package name: lisi)
check_pkg("kBET")
check_pkg("ggplot2")
check_pkg("dplyr")

# ----------------------------
# 1) Core metric function
# ----------------------------
score_embedding <- function(obj,
                            condition_col   = "condition",
                            celltype_col    = "final_population",
                            reduction       = "pca",
                            dims            = 1:30,
                            lisi_perplexity = 30,   # effective neighborhood size in LISI (not kNN k)
                            kbet_k0         = 15,   # neighborhood size for kBET
                            kbet_max_cells  = 5000, # subsample for speed (stratified by condition)
                            seed            = 1) {
  
  md <- obj@meta.data
  stopifnot(all(c(condition_col, celltype_col) %in% colnames(md)))
  stopifnot(reduction %in% names(obj@reductions))
  
  emb_all <- Seurat::Embeddings(obj, reduction = reduction)
  stopifnot(max(dims) <= ncol(emb_all))
  emb <- emb_all[, dims, drop = FALSE]
  
  # Keep only complete cases for the two labels
  keep <- stats::complete.cases(md[, c(condition_col, celltype_col), drop = FALSE])
  emb  <- emb[keep, , drop = FALSE]
  md   <- md[keep, , drop = FALSE]
  
  # Force factors (required for correct level counting)
  batch    <- factor(md[[condition_col]])
  celltype <- factor(md[[celltype_col]])
  
  B <- nlevels(batch)
  C <- nlevels(celltype)
  
  # ---- LISI (per-cell) ----
  # compute_lisi returns a data.frame with one column per label in label_colnames
  lisi_vals <- lisi::compute_lisi(
    X              = emb,
    meta_data      = md[, c(condition_col, celltype_col), drop = FALSE],
    label_colnames = c(condition_col, celltype_col),
    perplexity     = lisi_perplexity
  )
  
  iLISI <- lisi_vals[[condition_col]]
  cLISI <- lisi_vals[[celltype_col]]
  
  # Normalize iLISI to 0..1 where 0 = no mixing, 1 = maximal mixing given B batches
  iLISI_norm <- (iLISI - 1) / max(1, B - 1)
  
  # Convert cLISI to a "purity" score in 0..1 where 1 = perfectly pure neighborhoods
  # (cLISI = 1 means only one cell type locally; cLISI = C means full mixing)
  cPUR_norm <- (C - cLISI) / max(1, C - 1)
  
  stats_fun <- function(x) {
    c(
      mean   = mean(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE),
      p10    = unname(stats::quantile(x, 0.10, na.rm = TRUE)),
      p90    = unname(stats::quantile(x, 0.90, na.rm = TRUE))
    )
  }
  
  # ---- kBET ----
  set.seed(seed)
  idx <- seq_len(nrow(emb))
  if (length(idx) > kbet_max_cells) {
    by_batch <- split(idx, batch)
    sizes <- vapply(by_batch, length, integer(1))
    
    # Allocate approximately proportionally (guarantee at least 1 per batch if possible)
    alloc <- pmax(1L, floor(kbet_max_cells * sizes / sum(sizes)))
    # Trim if we overshot due to pmax(1)
    while (sum(alloc) > kbet_max_cells) {
      j <- which.max(alloc)
      if (alloc[j] > 1L) alloc[j] <- alloc[j] - 1L else break
    }
    
    idx <- sort(unlist(Map(function(ix, m) sample(ix, size = min(m, length(ix))), by_batch, alloc), use.names = FALSE))
  }
  
  X_sub     <- emb[idx, , drop = FALSE]
  batch_sub <- factor(batch[idx])
  
  # Use k0 and turn off heuristic so kBET does not silently change neighborhood size
  kbet_res <- tryCatch(
    kBET::kBET(
      df        = X_sub,
      batch     = batch_sub,
      k0        = kbet_k0,
      heuristic = FALSE,
      do.pca    = FALSE,
      plot      = FALSE,
      verbose   = FALSE
    ),
    error = function(e) e
  )
  
  kBET_mean_reject <- NA_real_
  kBET_p025 <- NA_real_
  kBET_p975 <- NA_real_
  
  if (!inherits(kbet_res, "error") && is.list(kbet_res)) {
    if (!is.null(kbet_res$stats) && "kBET.observed" %in% colnames(kbet_res$stats)) {
      obs <- kbet_res$stats[, "kBET.observed"]
      kBET_mean_reject <- mean(obs, na.rm = TRUE)
      kBET_p025 <- unname(stats::quantile(obs, 0.025, na.rm = TRUE))
      kBET_p975 <- unname(stats::quantile(obs, 0.975, na.rm = TRUE))
    } else if (!is.null(kbet_res$summary) && "kBET.observed" %in% names(kbet_res$summary)) {
      # Fallback: single observed rate from the summary
      kBET_mean_reject <- unname(kbet_res$summary[["kBET.observed"]])
    }
  }
  
  list(
    n_cells          = nrow(emb),
    n_batches        = B,
    n_celltypes      = C,
    iLISI_norm_stats = stats_fun(iLISI_norm),
    cPUR_norm_stats  = stats_fun(cPUR_norm),
    kBET_mean_reject = kBET_mean_reject,
    kBET_p025        = kBET_p025,
    kBET_p975        = kBET_p975
  )
}

# ----------------------------
# 2) Helpers
# ----------------------------
read_result_rdata <- function(path, object_name = NULL) {
  e <- new.env(parent = emptyenv())
  objs <- load(path, envir = e)
  
  if (!is.null(object_name)) {
    if (!object_name %in% objs) stop("Object '", object_name, "' not found in: ", path, call. = FALSE)
    return(e[[object_name]])
  }
  
  # Pick the first list that contains expected fields
  is_score <- function(x) {
    is.list(x) && all(c("iLISI_norm_stats", "cPUR_norm_stats", "kBET_mean_reject") %in% names(x))
  }
  hits <- objs[vapply(objs, function(nm) is_score(e[[nm]]), logical(1))]
  if (length(hits) == 0) stop("No score-like object found in: ", path, call. = FALSE)
  e[[hits[1]]]
}

# ----------------------------
# 3) Summarize and plot (three metrics)
# ----------------------------
summarize_scores <- function(res_list) {
  # res_list: named list of score_embedding() outputs, names are approach labels
  dplyr::bind_rows(lapply(names(res_list), function(nm) {
    res <- res_list[[nm]]
    data.frame(
      approach     = nm,
      n_cells      = res$n_cells,
      n_conditions = res$n_batches,
      n_celltypes  = res$n_celltypes,
      
      iLISI_med = unname(res$iLISI_norm_stats["median"]),
      iLISI_p10 = unname(res$iLISI_norm_stats["p10"]),
      iLISI_p90 = unname(res$iLISI_norm_stats["p90"]),
      
      cPUR_med  = unname(res$cPUR_norm_stats["median"]),
      cPUR_p10  = unname(res$cPUR_norm_stats["p10"]),
      cPUR_p90  = unname(res$cPUR_norm_stats["p90"]),
      
      kBET_mean = unname(res$kBET_mean_reject),
      kBET_p025 = unname(res$kBET_p025),
      kBET_p975 = unname(res$kBET_p975),
      
      stringsAsFactors = FALSE
    )
  }))
}

make_plot_df <- function(df_sum) {
  dplyr::bind_rows(
    dplyr::transmute(
      df_sum,
      approach,
      metric = "Condition mixing (iLISI, normalized)",
      center = iLISI_med, low = iLISI_p10, high = iLISI_p90
    ),
    dplyr::transmute(
      df_sum,
      approach,
      metric = "Cell-type purity (cLISI-derived, normalized)",
      center = cPUR_med, low = cPUR_p10, high = cPUR_p90
    ),
    dplyr::transmute(
      df_sum,
      approach,
      metric = "kBET rejection rate by condition",
      center = kBET_mean, low = kBET_p025, high = kBET_p975
    )
  )
}

plot_integration_metrics <- function(df_plot, approach_levels = NULL) {
  if (!is.null(approach_levels)) {
    df_plot$approach <- factor(df_plot$approach, levels = approach_levels)
  } else {
    df_plot$approach <- factor(df_plot$approach, levels = rev(unique(df_plot$approach)))
  }
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = approach, y = center)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = low, ymax = high),
      color = "grey60",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(color = "black", size = 2, na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ metric, scales = "fixed") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = NULL)
}


# res1 <- score_embedding(INTEGR_WEG0_PJ, condition_col = "condition", celltype_col = "CellType")
# res2 <- score_embedding(integrated_seurat_obj, condition_col = "condition", celltype_col = "possibly_final_anno")
# res3 <- score_embedding(result_obj, condition_col = "condition", celltype_col = "final_population")
# saveRDS(res1, "res1.rds"); saveRDS(res2, "res2.rds"); saveRDS(res3, "res3.rds")

# res_paths <- c(
#   "WT-only integrated (Harmony)"          = "G:/PhD_final/integration_presentation/res1.RData",
#   "All samples integrated (Harmony)"      = "G:/PhD_final/integration_presentation/res2.RData",
#   "All samples merged (no integration)"   = "G:/PhD_final/integration_presentation/res3.RData"
# )
# res_list <- lapply(res_paths, read_result_rdata)
# 
# df_sum  <- summarize_scores(res_list)
# df_plot <- make_plot_df(df_sum)
# 
# p <- plot_integration_metrics(df_plot, approach_levels = rev(names(res_paths)))
# 
# print(df_sum)
# print(p)