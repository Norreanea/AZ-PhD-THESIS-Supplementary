suppressPackageStartupMessages({
  library(Seurat)
  library(lisi)      # if unavailable on Windows, replace with the FNN-based fallback shown below
  library(kBET)
  library(dplyr)
  library(Matrix)
})



#load("D:/scRNA-seq/Gosia_obj/rep1/Schmed_integrated_SO_WEG0_FINAL_PJ.RData")
# An object of class Seurat 
# 35617 features across 22557 samples within 2 assays 
# Active assay: SCT (17808 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, integrated.harmony, umap.harmony
# head(INTEGR_WEG0_PJ$CellType)
# INTEGR_WEG0_PJ
# #pca, integrated.harmony, umap.harmony
# 
# p <- DimPlot(
#   INTEGR_WEG0_PJ,
#   group.by  = "CellType",
#   reduction = "umap.harmony",
#   label     = FALSE,
#   # we'll add labels in the next line
#   pt.size   = 1
#   #cols      = detailed_cols
# ) + NoLegend() + ggtitle("")
# 
# p <- LabelClusters(
#   plot     = p,
#   id       = "CellType",
#   repel    = TRUE,
#   box      = TRUE,
#   max.overlaps = Inf,
#   # <-- this is the key
#   size     = 3,
#   # label text size
#   label.size = 0.25           # box border thickness
# )
# 
# p
# cond <- sapply(strsplit(colnames(INTEGR_WEG0_PJ), "_"), `[`, 1)
# # Convert to uppercase for consistency
# cond <- toupper(cond)
# INTEGR_WEG0_PJ$condition <- cond
# unique(INTEGR_WEG0_PJ$condition)
# #INTEGR_WEG0_PJ$condition <- ifelse()
# p1 <- DimPlot(
#   INTEGR_WEG0_PJ,
#   group.by  = "condition",
#   reduction = "umap.harmony",
#   label     = FALSE,
#   # we'll add labels in the next line
#   pt.size   = 0.2
#   #cols      = detailed_cols
# ) + NoLegend() + ggtitle("")
# 
# p1 <- LabelClusters(
#   plot     = p1,
#   id       = "condition",
#   repel    = TRUE,
#   box      = TRUE,
#   max.overlaps = Inf,
#   # <-- this is the key
#   size     = 3,
#   # label text size
#   label.size = 0.25           # box border thickness
# )
# 
# p1
# 
# p2 <- DimPlot(
#   INTEGR_WEG0_PJ,
#   group.by  = "condition",
#   reduction = "pca",
#   label     = FALSE,
#   # we'll add labels in the next line
#   pt.size   = 0.2
#   #cols      = detailed_cols
# ) + NoLegend() + ggtitle("")
# 
# p2 <- LabelClusters(
#   plot     = p2,
#   id       = "condition",
#   repel    = TRUE,
#   box      = TRUE,
#   max.overlaps = Inf,
#   # <-- this is the key
#   size     = 3,
#   # label text size
#   label.size = 0.25           # box border thickness
# )
# 
# p2
# # 0) Common preprocessing (do once per approach, same HVGs and PCs)
# obj <- INTEGR_WEG0_PJ
# #obj <- NormalizeData(obj)
# obj <- FindVariableFeatures(obj, nfeatures = 3000)
# obj <- ScaleData(obj, features = VariableFeatures(obj))
# obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30)
# obj <- FindNeighbors(obj, dims = 1:30)
# IMPORTANT: do not compute metrics from UMAP

# 1) kBET (batch_key can be condition and, separately, timepoint)
#install_github('theislab/kBET') 

#emb <- Embeddings(obj, "pca")[,1:30]
#meta <- obj@meta.data
# emb  <- Embeddings(obj, "pca")[, 1:30, drop = FALSE]
# meta <- obj@meta.data
# meta$CellType <- factor(meta$CellType)
# 
# pc_r2 <- vapply(seq_len(ncol(emb)), function(j) {
#   fit <- lm(emb[, j] ~ meta$CellType)
#   stats::summary.lm(fit)$r.squared   # force the correct method
# }, numeric(1))
# 
# pcr_batch_mean <- mean(pc_r2, na.rm = TRUE)
# 
# pcr_batch_mean
# 
# set.seed(1)
# kbet_out <- kBET::kBET(emb, meta$condition, k0 = 15)  # lower rejection = better mixing
# 
# 
# library(lisi)
# library(Matrix)
# knn <- Seurat::Neighbors(obj)$nn  # or recompute with FNN::get.knn
# lisi_in <- list(X = emb, meta = meta[, c("condition","CellType")])
# lisi_vals <- compute_lisi(lisi_in$X, lisi_in$meta, c("condition","CellType"))  # columns iLISI, cLISI
# head(lisi_vals)
# B <- nlevels(factor(meta$condition))           # number of condition levels used
# C <- nlevels(factor(meta$CellType))    # number of cell-type levels
# 
# stats_fun <- function(x) c(mean=mean(x), median=median(x),
#                            p10=unname(quantile(x,0.10)), p90=unname(quantile(x,0.90)))
# 
# iLISI  <- lisi_vals[, "condition"]
# cLISI  <- lisi_vals[, "CellType"]
# 
# iLISI_norm <- (iLISI - 1) / (B - 1)   # 0 = no mixing, 1 = maximal mixing
# cLISI_norm <- (cLISI - 1) / (C - 1)   # 0 = perfect purity, 1 = maximal mixing
# 
# i_stats <- stats_fun(iLISI)      
# i_stats_norm <- stats_fun(iLISI_norm)
# c_stats <- stats_fun(cLISI)     
# c_stats_norm <- stats_fun(cLISI_norm)
# 
# i_stats
# i_stats_norm
# c_stats
# c_stats_norm


# Helper: compute LISI and kBET on a Seurat object
score_embedding <- function(obj,
                            condition_col = "condition",
                            celltype_col  = "final_population",
                            reduction     = "pca",
                            dims          = 1:30,
                            k_lisi        = 90,
                            k_kbet        = 15,
                            kbet_max_cells= 5000,
                            seed          = 1) {
  
  stopifnot(all(c(condition_col, celltype_col) %in% colnames(obj@meta.data)))
  emb  <- Seurat::Embeddings(obj, reduction)[, dims, drop = FALSE]
  meta <- obj@meta.data
  
  # --- LISI (per-cell) ---
  lisi_vals <- lisi::compute_lisi(emb,
                                  meta[, c(condition_col, celltype_col), drop = FALSE],
                                  c(condition_col, celltype_col))
  
  B <- nlevels(factor(meta[[condition_col]]))
  C <- nlevels(factor(meta[[celltype_col]]))
  iLISI      <- lisi_vals[, condition_col]
  cLISI      <- lisi_vals[, celltype_col]
  iLISI_norm <- (iLISI - 1) / max(1, B - 1)
  cLISI_norm <- (C - cLISI) / pmax(1, C - 1)
  
  # --- kBET (stratified subsample; positional call; robust arg names) ---
  set.seed(seed)
  idx <- seq_len(nrow(emb))
  if (length(idx) > kbet_max_cells) {
    by_cond <- split(idx, factor(meta[[condition_col]]))
    take <- lapply(by_cond, function(ix) {
      n <- length(ix); m <- ceiling(kbet_max_cells * n / length(idx))
      sample(ix, size = min(m, n))
    })
    idx <- sort(unlist(take, use.names = FALSE))
  }
  X     <- emb[idx, , drop = FALSE]
  batch <- factor(meta[[condition_col]][idx])
  
  kbet_formals <- names(formals(kBET::kBET))
  kbet_res <- tryCatch({
    if ("k0" %in% kbet_formals) {
      kBET::kBET(X, batch, k0 = k_kbet, plot = FALSE, do.pca = FALSE)
    } else if ("knn" %in% kbet_formals) {
      kBET::kBET(X, batch, knn = k_kbet, plot = FALSE, do.pca = FALSE)
    } else {
      kBET::kBET(X, batch, plot = FALSE, do.pca = FALSE)
    }
  }, error = function(e) e)
  
  kbet_mean_reject <- NA_real_
  kbet_detail <- NULL
  if (!inherits(kbet_res, "error")) {
    kbet_detail <- kbet_res$summary
    # compatible column picker
    cand <- c("kBET.observed", "observed", "rejection_rate")
    coln <- cand[cand %in% names(kbet_detail)][1]
    if (!is.null(coln)) kbet_mean_reject <- mean(kbet_detail[[coln]], na.rm = TRUE)
  }
  
  # --- summaries ---
  stats_fun <- function(x) c(mean = mean(x, na.rm = TRUE),
                             median = stats::median(x, na.rm = TRUE),
                             p10 = unname(stats::quantile(x, 0.10, na.rm = TRUE)),
                             p90 = unname(stats::quantile(x, 0.90, na.rm = TRUE)))
  
  out <- list(
    n_cells           = nrow(emb),
    n_batches         = B,
    n_celltypes       = C,
    iLISI_stats       = stats_fun(iLISI),
    iLISI_norm_stats  = stats_fun(iLISI_norm),
    cLISI_stats       = stats_fun(cLISI),
    cLISI_norm_stats  = stats_fun(cLISI_norm),
    kBET_mean_reject  = kbet_mean_reject,
    kBET_detail       = kbet_detail
  )
  return(out)
}

# Example across your three objects (adjust names/label columns to your data):
#res1 <- score_embedding(INTEGR_WEG0_PJ, condition_col = "condition", celltype_col = "CellType")
#res2 <- score_embedding(integrated_seurat_obj, condition_col = "condition", celltype_col = "possibly_final_anno")
# res1 <- score_embedding(INTEGR_WEG0_PJ,
#                         condition_col = "condition",
#                         celltype_col  = "CellType",
#                         reduction     = "pca",
#                         dims          = 1:30)
# save(res1,file="G:/PhD_final/integration_presentation/res1.RData")
# 
# length(unique(INTEGR_WEG0_PJ$CellType))
#res2 <- score_embedding(integrated_seurat_obj, condition_col = "condition", celltype_col = "possibly_final_anno")
#res3 <- score_embedding(result_obj, condition_col = "condition", celltype_col = "final_population")

# emb  <- Embeddings(obj, "pca")[, 1:30, drop=FALSE]
# X    <- model.matrix(~ factor(meta$CellType))   # fixed labels across all variants
# fit  <- lm.fit(X, emb)
# rss  <- colSums((emb - X %*% fit$coefficients)^2)
# centered <- sweep(emb, 2, colMeans(emb), "-")
# tss  <- colSums(centered^2)
# pc_r2 <- 1 - rss/tss
# 
# # weights proportional to PC variance (eigenvalues)
# w <- apply(centered, 2, var)     # or w <- obj@reductions$pca@stdev[1:30]^2
# r2_weighted_mean <- sum(w * pc_r2) / sum(w)
# r2_weighted_mean



load("G:/PhD_final/integration_presentation/res1.RData")
load("G:/PhD_final/integration_presentation/res2.RData")
load("G:/PhD_final/integration_presentation/res3.RData")
data.frame(
  approach = c("0h integrated", "All hpa integrated", "Merged"),
  iLISI_norm_stats = c(res1$iLISI_norm_stats, res2$iLISI_norm_stats, res3$iLISI_norm_stats),
  cLISI_norm_stats = c(res1$cLISI_norm_stats, res2$cLISI_norm_stats, res3$cLISI_norm_stats),
  kBET_mean_reject  = c(res1$kBET_mean_reject,  res2$kBET_mean_reject,  res3$kBET_mean_reject)
)


# iLISI_norm_stats["median"]: 0 = no local mixing by condition, 1 = maximal mixing given your number of conditions.
# cLISI_norm_stats["median"]: 0 = perfectly pure neighborhoods by cell type, larger values indicate mixing.
# kBET_mean_reject: average local rejection rate by condition; lower means better mixing by that factor.

# iLISI_norm_median: 0 means no local mixing by “condition”, 
# 1 means maximal mixing given your number of conditions.
# Values ≈ 0.3–0.5 suggest moderate mixing; values very close to 0 suggest strong separation by condition/timepoint.
# cLISI_norm_median: 0 means perfectly pure neighborhoods by cell type; 
# larger values indicate more local mixing of different types. You want this small.
# kBET_mean_reject: fraction of local tests rejecting the “well-mixed” null. 
# 0 means perfect mixing; ≥ 0.5 indicates strong separation by the tested factor. 
# You want this small when “condition” denotes an unwanted effect; 
# if “condition” includes true time biology, low rejection is not necessarily desirable.



# ---- build a clean summary table from res1/res2/res3 ----

pick_kbet_col <- function(kbet_detail) {
  if (is.null(kbet_detail)) return(NA_character_)
  cand <- c("kBET.observed", "observed", "rejection_rate")
  hit <- cand[cand %in% colnames(kbet_detail)][1]
  if (length(hit) == 0) NA_character_ else hit
}

summarize_one <- function(res, approach) {
  kcol <- pick_kbet_col(res$kBET_detail)
  
  k_p025 <- NA_real_
  k_p975 <- NA_real_
  if (!is.null(res$kBET_detail) && !is.na(kcol)) {
    rn <- rownames(res$kBET_detail)
    if ("2.5%" %in% rn)  k_p025 <- unname(res$kBET_detail["2.5%",  kcol])
    if ("97.5%" %in% rn) k_p975 <- unname(res$kBET_detail["97.5%", kcol])
  }
  
  data.frame(
    approach      = approach,
    n_cells       = res$n_cells,
    n_conditions  = res$n_batches,
    n_celltypes   = res$n_celltypes,
    
    iLISI_med     = unname(res$iLISI_norm_stats["median"]),
    iLISI_p10     = unname(res$iLISI_norm_stats["p10"]),
    iLISI_p90     = unname(res$iLISI_norm_stats["p90"]),
    
    # note: your cLISI_norm is a "purity" transform (higher = purer cell-type neighborhoods)
    cPUR_med      = unname(res$cLISI_norm_stats["median"]),
    cPUR_p10      = unname(res$cLISI_norm_stats["p10"]),
    cPUR_p90      = unname(res$cLISI_norm_stats["p90"]),
    
    kBET_mean     = unname(res$kBET_mean_reject),
    kBET_p025     = k_p025,
    kBET_p975     = k_p975,
    
    stringsAsFactors = FALSE
  )
}

df_sum <- dplyr::bind_rows(
  summarize_one(res1, "WT-only integrated (Harmony)"),
  summarize_one(res2, "All samples integrated (Harmony)"),
  summarize_one(res3, "Merged (no integration)")
)

# Optional: print a clean table for the thesis
df_sum


# ---- make a single, readable figure (three panels) ----

df_plot <- dplyr::bind_rows(
  df_sum |>
    dplyr::transmute(
      approach,
      metric = "Condition mixing (iLISI, normalized)",
      center = iLISI_med,
      low    = iLISI_p10,
      high   = iLISI_p90
    ),
  df_sum |>
    dplyr::transmute(
      approach,
      metric = "Cell-type purity (cLISI-derived, normalized)",
      center = cPUR_med,
      low    = cPUR_p10,
      high   = cPUR_p90
    ),
  df_sum |>
    dplyr::transmute(
      approach,
      metric = "kBET rejection rate by condition",
      center = kBET_mean,
      low    = kBET_p025,
      high   = kBET_p975
    )
)

df_plot$approach <- factor(df_plot$approach, levels = df_sum$approach)
#?facet_wrap
#?aes
# ggplot2::ggplot(df_plot, ggplot2::aes(x = approach, y = center, color ="black")) +
#   ggplot2::geom_pointrange(ggplot2::aes(ymin = low, ymax = high, color="grey"), na.rm = TRUE) +
#   ggplot2::coord_flip() +
#   ggplot2::facet_wrap(~ metric, scales = "fixed") +
#   ggplot2::theme_classic() +
#   ggplot2::labs(x = NULL, y = NULL)
df_plot$approach <- as.character(df_plot$approach)
df_plot$approach <- ifelse(df_plot$approach=="Merged (no integration)","All samples merged (no integration)",df_plot$approach)
df_plot$approach <- factor(df_plot$approach, levels = c("All samples merged (no integration)","All samples integrated (Harmony)","WT-only integrated (Harmony)"))
df_plot$approach <- factor(df_plot$approach, levels = rev(unique(df_plot$approach)))

ggplot2::ggplot(df_plot, ggplot2::aes(x = approach, y = center)) +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = low, ymax = high),
                           color = "grey60", na.rm = TRUE) +
  ggplot2::geom_point(color = "black", size = 2, na.rm = TRUE) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ metric, scales = "fixed") +
  ggplot2::theme_classic() +
  ggplot2::labs(x = NULL, y = NULL)



df_heat <- df_plot |>
  dplyr::mutate(
    metric = factor(metric, levels = unique(metric)),
    lbl = dplyr::if_else(
      is.na(low) | is.na(high),
      sprintf("%.2f", center),
      sprintf("%.2f\n[%.2f–%.2f]", center, low, high)
    ),
    # Color score: make "green = better" consistent across metrics
    fill_score = dplyr::if_else(
      metric == "kBET rejection rate by condition",
      1 - center,   # lower kBET = better -> higher score
      center        # higher iLISI/cPUR = better
    )
  )
df_heat$metric <- as.character(df_heat$metric)
df_heat$metric <- ifelse(df_heat$metric=="Condition mixing (iLISI, normalized)","Condition mixing",df_heat$metric)
df_heat$metric <- ifelse(df_heat$metric=="Cell-type purity (cLISI-derived, normalized)","Cell-type purity",df_heat$metric)
df_heat$metric <- ifelse(df_heat$metric=="kBET rejection rate by condition","kBET rejection rate",df_heat$metric)
p_heat <- ggplot2::ggplot(df_heat, ggplot2::aes(x = metric, y = approach, fill = fill_score)) +
  ggplot2::geom_tile() +
  ggplot2::geom_label(
    ggplot2::aes(label = lbl),
    size = 4,
    color = "black",
    fill = "white",
    label.size = 0.3,
    label.padding = grid::unit(0.12, "lines")
  ) +
  ggplot2::scale_fill_gradient(
    low = "pink",
    high = "skyblue",
    limits = c(0, 1),
    name = "Score"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(vjust = 0.5),
    legend.position = "right",
    text = element_text(size = 15)
  )

p_heat
# If you want to save:
# ggplot2::ggsave("integration_metrics.png", p, width = 10, height = 5, dpi = 300)
