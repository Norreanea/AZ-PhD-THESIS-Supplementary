# Auto-annotation pipeline 

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(openxlsx)
  library(igraph)
})
suppressPackageStartupMessages(library(future))
plan(sequential)

if (.Platform$OS.type == "windows" && exists("memory.limit")) {
  try(suppressWarnings(memory.limit(size = 56000)), silent = TRUE)
}

# ---------------------------
# Config (tuned defaults)
# ---------------------------
cfg <- list(
  paths = list(
    seurat_rdata   = "D:/scRNA-seq/AZ_final_obj/seurat_obj_new.RData",
    matrix_rds     = NULL,
    markers_xlsx   = "G:/PhD_final/cell_markers_curated_new_new.xlsx",
    anno_rdata     = "E:/Stringtie_anno/SM_anno/final/final_final/pfam_swiss_ncbi_merged_only_genes_dedup.RData",
    out_xlsx       = sprintf(
      "G:/PhD_final/auto_annotation_%s.xlsx",
      format(Sys.time(), "%Y%m%d_%H%M")
    )
  ),
  ckpt_dir        = "G:/PhD_final/sncRNA/.auto_annot_ckpts",
  use_qs          = TRUE,
  qs_preset       = "balanced",
  
  base_assay      = "SCT",
  pca_name        = "pca.auto",
  umap_name       = "umap.auto",
  max_pcs         = 60L,
  variance_cut    = 0.90,
  knee_smooth     = 5L,
  
  target_n_clusters = 60L,
  k_grid            = c(5L, 8L, 10L, 15L, 20L),
  res_grid          = c(seq(0.4, 2.0, by = 0.2), 1.8, 2.0),
  res_init          = 0.6,
  res_max           = 10,
  grid_max_steps    = 20,
  
  tiny_frac_cut     = 0.015,
  agree_cut         = 0.75,
  
  sub_npcs          = 30L,
  seed              = 42L,
  sub_k_grid        = c(10L, 15L, 20L, 25L),
  sub_res_grid      = seq(0.5, 3.5, by = 0.25),
  sub_min_cells_for_split = 30L,
  sub_max_children  = 6L,
  sub_min_child_n    = 10L,
  sub_min_child_prop = 0.005,
  
  annot_features_max    = 1500L,
  annot_skip_deg        = TRUE,
  deg_subsample_per_ident = 2000L,
  
  enable_ucell     = TRUE,
  ucell_min_genes  = 3L,
  ucell_cells_per_cluster = 1000L,
  ucell_max_signatures   = 300L,
  ucell_ncores     = 1L,
  
  write_round1_degs = TRUE
)
set.seed(cfg$seed)

# Optional packages
has_ucell     <- requireNamespace("UCell", quietly = TRUE)
has_fgsea     <- requireNamespace("fgsea", quietly = TRUE)
has_cellmanam <- requireNamespace("CellMaNam", quietly = TRUE)
has_qs        <- requireNamespace("qs", quietly = TRUE)
has_digest    <- requireNamespace("digest", quietly = TRUE)

# ---------------------------
# Helpers
# ---------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
trim <- function(x)
  gsub("^\\s+|\\s+$", "", x)
canon_cluster <- function(v) {
  v <- as.character(v)
  v <- trimws(v)
  v <- sub("^X([0-9]+)$", "g\\1", v)
  v <- sub("^([0-9]+)$", "g\\1", v)
  v
}
sanitize_sheet <- function(x) {
  x <- gsub("[\\*\\?/\\\\\\[\\]:]", "_", x)
  x <- substr(x, 1, 31)
  make.unique(x)
}
pick_layer_arg <- function() {
  if ("layer" %in% names(formals(Seurat::FindAllMarkers)))
    "layer"
  else
    "slot"
}
apply_mapping <- function(keys, map_named) {
  out <- rep(NA_character_, length(keys))
  m <- match(keys, names(map_named))
  hit <- !is.na(m)
  out[hit] <- unname(map_named[m[hit]])
  out
}
as_chr_collapse <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- unique(na.omit(as.character(x)))
  if (!length(x)) "" else paste(x, collapse = "; ")
}


# IO/ckpt ------------------------------------------------------------
CKPT_DIR <- cfg$ckpt_dir
.dir_ok <- function() {
  dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
  TRUE
}
ckpt_path_qs  <- function(stage)
  file.path(CKPT_DIR, paste0("auto_annot_ckpt_", stage, ".qs"))
ckpt_path_rds <- function(stage)
  file.path(CKPT_DIR, paste0("auto_annot_ckpt_", stage, ".rds"))
ckpt_has  <- function(stage)
  file.exists(ckpt_path_qs(stage)) ||
  file.exists(ckpt_path_rds(stage))
ckpt_save <- function(stage, value) {
  .dir_ok()
  if (has_qs &&
      isTRUE(cfg$use_qs))
    qs::qsave(value, ckpt_path_qs(stage), preset = cfg$qs_preset)
  else
    saveRDS(value, ckpt_path_rds(stage))
}
ckpt_load <- function(stage) {
  if (file.exists(ckpt_path_qs(stage)))
    return(qs::qread(ckpt_path_qs(stage)))
  readRDS(ckpt_path_rds(stage))
}
wb_load_or_new <- function(path) if (file.exists(path)) openxlsx::loadWorkbook(path) else openxlsx::createWorkbook()
wb_save <- function(wb, path)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
ckpt_update <- function(stage, st, ...) {
  up <- list(...)
  for (nm in names(up))
    st[[nm]] <- up[[nm]]
  ckpt_save(stage, st)
}

# Memory diet ------------------------------------------------------------
diet_for_checkpoint <- function(obj,
                                keep_assay  = cfg$base_assay,
                                keep_reduc  = c(cfg$pca_name, cfg$umap_name, paste0(cfg$umap_name, ".v2")),
                                keep_graphs = character(),
                                drop_counts = TRUE,
                                drop_scale  = TRUE) {
  DefaultAssay(obj) <- keep_assay
  obj@assays     <- obj@assays[intersect(names(obj@assays), keep_assay)]
  obj@reductions <- obj@reductions[intersect(names(obj@reductions), keep_reduc)]
  obj@graphs     <- obj@graphs[intersect(names(obj@graphs), keep_graphs)]
  
  v <- tryCatch(
    utils::packageVersion("SeuratObject"),
    error = function(e)
      package_version("4.0.0")
  )
  if (v >= package_version("5.0.0")) {
    Seurat::DietSeurat(
      obj,
      assays = names(obj@assays),
      dimreducs = names(obj@reductions),
      graphs = names(obj@graphs),
      layers = setNames(list("data"), names(obj@assays))
    )
  } else {
    Seurat::DietSeurat(
      obj,
      assays = names(obj@assays),
      counts = !drop_counts,
      data = TRUE,
      scale.data = !drop_scale,
      dimreducs = names(obj@reductions),
      graphs = names(obj@graphs),
      features = NULL
    )
  }
}

# Stats helpers ----------------------------------------------------------
.named_or_fail <- function(x, expect_names, method_name) {
  if (!length(x))
    return(x)
  if (is.null(names(x)) || any(is.na(names(x)) | names(x) == "")) {
    if (!missing(expect_names) && length(x) == length(expect_names)) {
      names(x) <- as.character(expect_names)
    } else {
      stop(
        sprintf(
          "Annotation method produced an unnamed vector; method=%s len=%d",
          method_name,
          length(x)
        )
      )
    }
  }
  names(x) <- canon_cluster(names(x))
  x
}
.run_annot_method <- function(label,
                              fun,
                              expect_clusters,
                              ckpt_stage_tag = "annot") {
  message(sprintf("[annot] %s: running...", label))
  v <- tryCatch(
    fun(),
    error = function(e) {
      message(sprintf("[annot] %s: ERROR -> %s", label, conditionMessage(e)))
      return(setNames(character(0), character(0)))
    }
  )
  v <- .named_or_fail(v, expect_clusters, label)
  if (length(v))
    v <- v[names(v) %in% as.character(expect_clusters)]
  tag <- paste0(ckpt_stage_tag, "_", label)
  try(ckpt_save(
    tag,
    list(
      method = label,
      result = v,
      expect_clusters = as.character(expect_clusters),
      saved_at = Sys.time()
    )
  ), silent = TRUE)
  message(sprintf("[annot] %s: %d label(s).", label, length(v)))
  v
}

choose_pcs_by_knee <- function(stdev,
                               max_pcs = 50L,
                               variance_cut = 0.90,
                               smooth_k = 5L,
                               min_pcs = 25L) {
  stdev <- stdev[is.finite(stdev) & stdev > 0]
  stdev <- stdev[seq_len(min(length(stdev), max_pcs))]
  if (!length(stdev))
    return(min_pcs)
  
  var_ratio <- stdev^2 / sum(stdev^2)
  cumvar <- cumsum(var_ratio)
  ceil <- which(cumvar >= variance_cut)[1]
  if (is.na(ceil))
    ceil <- length(var_ratio)
  
  y <- var_ratio
  if (length(y) >= (smooth_k * 2 + 1)) {
    y <- stats::filter(y, rep(1 / (smooth_k * 2 + 1), smooth_k * 2 + 1), sides = 2)
    y[is.na(y)] <- var_ratio[is.na(y)]
  }
  d2 <- diff(y, differences = 2)
  knee <- which.min(d2) + 1L
  
  pcs <- max(min_pcs, min(max(10L, knee), ceil))
  pcs <- min(pcs, length(stdev))
  pcs
}

avg_by_group_sparse <- function(X, groups) {
  stopifnot(inherits(X, "dgCMatrix"))
  groups <- droplevels(factor(groups))
  G <- Matrix::sparse.model.matrix(~ groups - 1)
  sums <- X %*% G
  n_per <- Matrix::colSums(G)
  Dinv <- Matrix::Diagonal(x = as.numeric(1 / n_per))
  avg <- sums %*% Dinv
  colnames(avg) <- levels(groups)
  avg
}
get_graph_modularity <- function(obj, graph.name, membership) {
  S <- obj@graphs[[graph.name]]
  if (is.null(S))
    return(NA_real_)
  if (!methods::is(S, "dgCMatrix"))
    S <- as(S, "dgCMatrix")
  S <- Matrix::drop0((S + Matrix::t(S)) / 2)
  Matrix::diag(S) <- 0
  trip <- as.data.frame(Matrix::summary(S))
  trip <- trip[trip$i < trip$j &
                 trip$x > 0, c("i", "j", "x"), drop = FALSE]
  if (!nrow(trip))
    return(NA_real_)
  vnames <- colnames(S)
  if (is.null(vnames))
    vnames <- seq_len(ncol(S))
  g <- igraph::graph_from_data_frame(
    data.frame(
      from = vnames[trip$i],
      to = vnames[trip$j],
      weight = trip$x
    ),
    directed = FALSE,
    vertices = data.frame(name = vnames)
  )
  mvec <- membership
  if (!is.null(names(mvec)))
    mvec <- mvec[V(g)$name]
  else
    names(mvec) <- V(g)$name
  memb <- as.integer(factor(mvec))
  igraph::modularity(g, membership = memb, weights = igraph::E(g)$weight)
}
score_solution <- function(obj,
                           graph.name,
                           idents,
                           target_n,
                           tiny_cut) {
  memb <- as.character(idents)
  nclu <- length(unique(memb))
  n    <- length(memb)
  tab  <- sort(table(memb), decreasing = TRUE)
  tiny_frac <- if (length(tab))
    sum(tab < (tiny_cut * n)) / length(tab)
  else
    1
  mod <- suppressWarnings(get_graph_modularity(obj, graph.name, memb))
  if (!is.finite(mod) || mod < 1e-6) {
    return(data.frame(
      n_clusters = nclu,
      modularity = mod,
      tiny_frac = tiny_frac,
      score = -Inf
    ))
  }
  if (is.na(target_n))
    close <- 0
  else {
    close <- 1 - (abs(nclu - target_n) / max(target_n, nclu))
    close <- max(0, min(close, 1))
  }
  tiny_pen <- pmin(0.4, tiny_frac * 0.8)
  score <- 0.70 * close + 0.25 * mod - 0.05 * tiny_pen
  data.frame(
    n_clusters = nclu,
    modularity = mod,
    tiny_frac = tiny_frac,
    score = score
  )
}

# Marker prep ------------------------------------------------------------
prepare_marker_ref <- function(cell_markers_df,
                               gene_col = "Markers_positive_SMESG",
                               general_col = "Cell_population_general",
                               detailed_col = "Cell_population_detailed",
                               neg_col = NULL,
                               weight_col = NULL) {
  gene_col <- match.arg(gene_col)
  df <- as.data.frame(cell_markers_df, stringsAsFactors = FALSE)
  norm <- function(v)
    trim(as.character(v))
  gene     <- norm(df[[gene_col]])
  general  <- norm(df[[general_col]])
  detailed <- norm(df[[detailed_col]])
  neg_vec <- if (!is.null(neg_col) &&
                 neg_col %in% names(df))
    as.logical(df[[neg_col]])
  else
    FALSE
  weight_vec <- if (!is.null(weight_col) &&
                    weight_col %in% names(df))
    suppressWarnings(as.numeric(df[[weight_col]]))
  else
    NA_real_
  general[general %in% c("Protonephridia", "Protonephridia ")] <- "Protonephridia"
  detailed[detailed %in% c("Protonephridial tubule precursor ",
                           "Protonephridia tubule precursor")] <- "Protonephridial tubule precursor"
  pos <- tibble::tibble(
    gene = gene,
    general = general,
    detailed = detailed,
    neg = neg_vec,
    weight = weight_vec
  ) %>%
    dplyr::filter(gene != "" &
                    general != "" &
                    detailed != "") %>% dplyr::distinct()
  ref_general  <- pos %>% dplyr::transmute(gene, final_cluster = general, weight, neg)
  ref_detailed <- pos %>% dplyr::transmute(gene,
                                           final_cluster = detailed,
                                           parent_general = general,
                                           weight,
                                           neg)
  list(general = ref_general, detailed = ref_detailed)
}
.bg_genes <- function(obj)
  rownames(Seurat::GetAssayData(obj, assay = cfg$base_assay, layer = "data"))

# DEG robust + cache ------------------------------------------------------
compute_degs_robust <- function(obj, group_col, features_whitelist = NULL) {
  DefaultAssay(obj) <- cfg$base_assay
  stopifnot(group_col %in% colnames(obj@meta.data))
  pick <- function(df, cands, default) {
    for (nm in cands)
      if (nm %in% names(df))
        return(df[[nm]])
    if (is.function(default))
      return(default())
    rep(default, nrow(df))
  }
  if (is.null(features_whitelist)) {
    bg   <- rownames(Seurat::GetAssayData(obj, assay = cfg$base_assay, layer = "data"))
    hvgs <- tryCatch(
      VariableFeatures(obj),
      error = function(e)
        character(0)
    )
    if (!length(hvgs)) {
      obj <- FindVariableFeatures(
        obj,
        assay = cfg$base_assay,
        nfeatures = 3000,
        verbose = FALSE
      )
      hvgs <- VariableFeatures(obj)
    }
    misc_markers <- tryCatch(
      unique(unlist(obj@misc$marker_genes)),
      error = function(e)
        character(0)
    )
    features_whitelist <- unique(intersect(union(hvgs, misc_markers), bg))
    if (!length(features_whitelist))
      features_whitelist <- hvgs[hvgs %in% bg]
    if (!length(features_whitelist))
      features_whitelist <- bg
  }
  layer_or_slot <- pick_layer_arg()
  old_id <- Idents(obj)
  on.exit(Idents(obj) <- old_id, add = TRUE)
  Idents(obj) <- obj[[group_col]][, 1]
  lv <- levels(Idents(obj))
  if (!length(lv)) {
    message("compute_degs_robust: no levels in '", group_col, "'.")
    return(
      tibble::tibble(
        cluster = character(),
        gene = character(),
        avg_log2FC = double(),
        p_val_adj = double()
      )
    )
  }
  pieces <- lapply(lv, function(cl) {
    args <- list(
      object = obj,
      ident.1 = cl,
      only.pos = TRUE,
      min.pct = 0.20,
      logfc.threshold = 0.25,
      test.use = "wilcox",
      verbose = FALSE,
      random.seed = cfg$seed,
      features = features_whitelist
    )
    if (!is.null(cfg$deg_subsample_per_ident))
      args$max.cells.per.ident <- cfg$deg_subsample_per_ident
    if (layer_or_slot == "layer")
      args$layer <- "data"
    else
      args$slot <- "data"
    fm <- tryCatch(
      do.call(Seurat::FindMarkers, args),
      error = function(e)
        NULL
    )
    if (is.null(fm) || !nrow(fm))
      return(NULL)
    tibble::tibble(
      cluster    = canon_cluster(cl),
      gene       = rownames(fm),
      avg_log2FC = suppressWarnings(pick(
        fm, c("avg_log2FC", "avg_logFC", "log2FC"), NA_real_
      )),
      p_val_adj  = suppressWarnings(pick(
        fm,
        c("p_val_adj", "p_val.adj", "p_val_adj_fdr", "p_val"),
        1
      ))
    ) %>% dplyr::filter(is.finite(avg_log2FC),
                        is.finite(p_val_adj),
                        p_val_adj < 0.05,
                        avg_log2FC > 0)
  })
  res <- dplyr::bind_rows(pieces)
  if (!nrow(res)) {
    message("compute_degs_robust: no DEGs passed filters.")
    return(
      tibble::tibble(
        cluster = character(),
        gene = character(),
        avg_log2FC = double(),
        p_val_adj = double()
      )
    )
  }
  res
}
.safe_hash <- function(x) {
  x <- paste(x, collapse = "|")
  if (has_digest)
    digest::digest(x, algo = "xxhash64")
  else
    sprintf("h%08x", abs(as.integer(sum(utf8ToInt(
      x
    ))) %% 2^31))
}
.deg_ckpt_tag <- function(obj, group_col, features = NULL) {
  memb <- canon_cluster(as.character(obj[[group_col]][, 1]))
  sz   <- sort(as.integer(table(memb)), decreasing = TRUE)
  features <- sort(unique(as.character(features %||% character(0))))
  feat_stub <- features[seq_len(min(200L, length(features)))]
  .safe_hash(c(
    sprintf("n=%d", length(memb)),
    sprintf("k=%d", length(sz)),
    paste0("sz:", paste(sz, collapse = ",")),
    sprintf("f=%d", length(features)),
    paste0("feat:", paste(feat_stub, collapse = ","))
  ))
}
.deg_ckpt_file <- function(prefix, tag)
  file.path(CKPT_DIR, sprintf("deg_%s_%s.rds", prefix, tag))
deg_ckpt_save <- function(prefix,
                          obj,
                          group_col,
                          degs,
                          features = NULL,
                          extra = list()) {
  .dir_ok()
  tag <- .deg_ckpt_tag(obj, group_col, features)
  saveRDS(c(
    list(
      saved_at = Sys.time(),
      prefix = prefix,
      group_col = group_col,
      tag = tag,
      features = features,
      degs = degs
    ),
    extra
  ), .deg_ckpt_file(prefix, tag))
}
deg_ckpt_load <- function(prefix, obj, group_col, features = NULL) {
  tag <- .deg_ckpt_tag(obj, group_col, features)
  f <- .deg_ckpt_file(prefix, tag)
  if (!file.exists(f))
    return(NULL)
  x <- readRDS(f)
  if (!is.list(x) || is.null(x$degs))
    return(NULL)
  x$degs
}

# Neighbors/Clustering ----------------------------------------------------
run_neighbors_if_needed <- function(obj, dims, k, reduction = cfg$pca_name) {
  gname <- paste0(cfg$base_assay, "_snn_k", k)
  if (is.null(obj@graphs[[gname]])) {
    obj <- FindNeighbors(
      obj,
      reduction = reduction,
      dims = 1:dims,
      k.param = k,
      graph.name = gname,
      verbose = FALSE
    )
  }
  obj
}
safe_findclusters <- function(obj, graph.name, resolution, seed = cfg$seed) {
  ok <- FALSE
  res <- NULL
  for (alg in c(4, 3, 1)) {
    res <- try(FindClusters(
      obj,
      graph.name = graph.name,
      resolution = resolution,
      algorithm = alg,
      random.seed = seed,
      verbose = FALSE
    ),
    silent = TRUE)
    if (!inherits(res, "try-error")) {
      ok <- TRUE
      break
    }
  }
  if (!ok)
    stop("FindClusters failed for graph=",
         graph.name,
         " res=",
         resolution)
  res
}

# Knee-based resolution picker -------------------------------------------
pick_res_by_knee <- function(res_vals, nclu_vals, target_n = NA_integer_) {
  stopifnot(length(res_vals) == length(nclu_vals), length(res_vals) > 0)
  o <- order(res_vals)
  x <- as.numeric(res_vals[o])
  y <- as.numeric(nclu_vals[o])
  y <- cummax(y)
  dy  <- c(NA_real_, diff(y))
  d2  <- c(NA_real_, diff(dy))
  elbow_ix <- suppressWarnings(which.max(replace(-d2, is.na(d2), -Inf)))
  if (length(elbow_ix) == 0 ||
      is.infinite(elbow_ix) || is.na(elbow_ix) || elbow_ix < 1) {
    if (is.na(target_n))
      return(x[ceiling(length(x) / 2)])
    return(x[which.min(abs(y - target_n))])
  }
  x[elbow_ix]
}

# Grid search with stability & memory care ---------------------------
grid_search_clusters <- function(obj,
                                 dims,
                                 target_n = cfg$target_n_clusters,
                                 k_grid = c(5L, 8L, 10L, 15L, 20L),
                                 res_init = 0.6,
                                 res_max  = 10,
                                 max_steps = 20) {
  attempts <- list()
  best <- NULL
  best_membership <- NULL
  best_k <- NULL
  best_r <- NULL
  
  for (k in k_grid) {
    obj_k <- run_neighbors_if_needed(obj, dims, k)
    gname <- paste0(cfg$base_assay, "_snn_k", k)
    
    eval_r <- function(r) {
      x <- safe_findclusters(
        obj_k,
        graph.name = gname,
        resolution = r,
        seed = cfg$seed
      )
      memb <- as.character(Idents(x))
      names(memb) <- colnames(x)
      rec  <- score_solution(x, gname, memb, target_n, cfg$tiny_frac_cut)
      rec$k <- k
      rec$resolution <- r
      list(rec = rec,
           memb = memb,
           obj = x)
    }
    
    r_lo <- res_init
    e_lo <- eval_r(r_lo)
    n_lo <- e_lo$rec$n_clusters
    attempts[[length(attempts) + 1]] <- e_lo$rec
    
    if (is.na(target_n)) {
      for (r in seq(res_init, min(res_max, res_init + 4), by = 0.4)) {
        e <- eval_r(r)
        attempts[[length(attempts) + 1]] <- e$rec
        if (is.null(best) ||
            e$rec$score > best$score) {
          best <- e$rec
          best_membership <- e$memb
          best_k <- k
          best_r <- r
        }
      }
      next
    }
    
    r_hi <- r_lo
    e_hi <- e_lo
    n_hi <- n_lo
    step <- 0
    while (n_hi < target_n && r_hi < res_max && step < max_steps) {
      r_hi <- r_hi * 1.5
      if (r_hi <= r_lo)
        r_hi <- r_lo + 0.2
      e_hi <- eval_r(r_hi)
      attempts[[length(attempts) + 1]] <- e_hi$rec
      n_hi <- e_hi$rec$n_clusters
      step <- step + 1
    }
    
    cand_list <- list(e_lo, e_hi)
    if (n_hi < target_n) {
      e_closest <- cand_list[[which.min(abs(c(n_lo, n_hi) - target_n))]]
      if (is.null(best) || e_closest$rec$score > best$score) {
        best <- e_closest$rec
        best_membership <- e_closest$memb
        best_k <- k
        best_r <- e_closest$rec$resolution
      }
      next
    }
    
    l_r <- r_lo
    l_e <- e_lo
    l_n <- n_lo
    h_r <- r_hi
    h_e <- e_hi
    h_n <- n_hi
    
    it <- 0
    while (it < max_steps && (abs(h_r - l_r) > 0.05)) {
      it <- it + 1
      m_r <- (l_r + h_r) / 2
      m_e <- eval_r(m_r)
      attempts[[length(attempts) + 1]] <- m_e$rec
      m_n <- m_e$rec$n_clusters
      if (m_n < target_n) {
        l_r <- m_r
        l_e <- m_e
        l_n <- m_n
      } else {
        h_r <- m_r
        h_e <- m_e
        h_n <- m_n
      }
    }
    
    final_e <- if (abs(l_n - target_n) <= abs(h_n - target_n))
      l_e
    else
      h_e
    if (is.null(best) || final_e$rec$score > best$score) {
      best <- final_e$rec
      best_membership <- final_e$memb
      best_k <- k
      best_r <- final_e$rec$resolution
    }
    gc(FALSE)
  }
  
  stopifnot(!is.null(best))
  cl_fac <- factor(best_membership, levels = unique(best_membership))
  names(cl_fac) <- names(best_membership)
  Idents(obj) <- cl_fac[colnames(obj)]
  obj$seurat_clusters <- Idents(obj)
  
  keep_graph <- paste0(cfg$base_assay, "_snn_k", best_k)
  obj@graphs <- obj@graphs[intersect(names(obj@graphs), keep_graph)]
  
  diag <- dplyr::bind_rows(attempts)
  obj@misc$grid_diag <- diag
  
  list(obj = obj,
       stats = best,
       keep_graph = keep_graph)
}

# Annotation methods ------------------------------------------------------
annotate_avgexp_matrix <- function(avg_mat,
                                   marker_ref,
                                   curated_genes,
                                   top_n = 10) {
  stopifnot(is.matrix(avg_mat) || inherits(avg_mat, "Matrix"))
  feats <- intersect(curated_genes, rownames(avg_mat))
  if (!length(feats))
    return(setNames(character(0), character(0)))
  Z <- t(scale(t(as.matrix(avg_mat[feats, , drop = FALSE]))))
  Z[is.na(Z)] <- 0
  res <- vapply(seq_len(ncol(Z)), function(i) {
    top <- head(names(sort(Z[, i], decreasing = TRUE)), top_n)
    tab <- dplyr::filter(marker_ref, gene %in% top) %>% dplyr::count(final_cluster, sort = TRUE)
    if (nrow(tab) == 0)
      "Unknown"
    else
      tab$final_cluster[1]
  }, FUN.VALUE = character(1))
  names(res) <- colnames(Z)
  res
}
annotate_hypergeom <- function(seurat_clusters_DEG,
                               marker_ref,
                               bg_genes) {
  bg_genes <- unique(bg_genes)
  if (!nrow(seurat_clusters_DEG))
    return(setNames(character(0), character(0)))
  clusters <- unique(seurat_clusters_DEG$cluster)
  res <- sapply(clusters, function(cl) {
    cl_genes <- intersect(unique(seurat_clusters_DEG$gene[seurat_clusters_DEG$cluster == cl]), bg_genes)
    total <- length(bg_genes)
    df <- do.call(rbind, lapply(unique(marker_ref$final_cluster), function(ct) {
      ct_genes <- intersect(unique(marker_ref$gene[marker_ref$final_cluster == ct]), bg_genes)
      overlap <- length(intersect(cl_genes, ct_genes))
      m <- length(ct_genes)
      n <- total - m
      k <- length(cl_genes)
      pval <- stats::phyper(overlap - 1, m, n, k, lower.tail = FALSE)
      data.frame(cell_type = ct, pval = pval)
    }))
    df$p_adj <- p.adjust(df$pval, method = "BH")
    df$cell_type[which.min(df$p_adj)]
  })
  names(res) <- clusters
  res
}
annotate_majority <- function(seurat_clusters_DEG, marker_ref) {
  if (!nrow(seurat_clusters_DEG))
    return(setNames(character(0), character(0)))
  df <- dplyr::inner_join(seurat_clusters_DEG,
                          marker_ref,
                          by = "gene",
                          relationship = "many-to-many") %>%
    dplyr::group_by(cluster, final_cluster) %>% dplyr::summarise(n = dplyr::n(), .groups =
                                                                   "drop") %>%
    dplyr::group_by(cluster) %>% dplyr::slice_max(n, n = 1, with_ties =
                                                    FALSE)
  setNames(df$final_cluster, df$cluster)
}
annotate_logfc <- function(seurat_clusters_DEG,
                           marker_ref,
                           gene_specificity = NULL,
                           use_padj_weight = TRUE,
                           neg_penalty = 0.2) {
  if (!nrow(seurat_clusters_DEG))
    return(setNames(character(0), character(0)))
  df <- dplyr::inner_join(seurat_clusters_DEG,
                          marker_ref,
                          by = "gene",
                          relationship = "many-to-many")
  if (use_padj_weight &&
      "p_val_adj" %in% names(seurat_clusters_DEG)) {
    df <- dplyr::mutate(df, w = pmax(0, -log10(pmin(
      p_val_adj, 1e-300
    ))))
  } else
    df$w <- 1
  if (!is.null(gene_specificity))
    df$spec <- pmax(0.2, gene_specificity[match(df$gene, names(gene_specificity))])
  else
    df$spec <- 1
  if ("neg" %in% names(df) &&
      any(!is.na(df$neg)))
    df$neg_w <- ifelse(isTRUE(df$neg), -neg_penalty, 0)
  else
    df$neg_w <- 0
  df <- df %>%
    dplyr::group_by(cluster, final_cluster) %>%
    dplyr::summarise(
      score = stats::weighted.mean(pmax(0, avg_log2FC) * spec + neg_w, w, na.rm =
                                     TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(cluster) %>% dplyr::slice_max(score, n = 1, with_ties =
                                                    FALSE)
  setNames(df$final_cluster, df$cluster)
}
annotate_cellmanam <- function(seurat_obj,
                               seurat_clusters_DEG,
                               marker_ref,
                               top_n = 2,
                               p_val = 0.01,
                               level = 1) {
  if (!isTRUE(has_cellmanam))
    return(setNames(character(0), character(0)))
  if (is.null(seurat_clusters_DEG) ||
      !nrow(seurat_clusters_DEG))
    return(setNames(character(0), character(0)))
  DefaultAssay(seurat_obj) <- cfg$base_assay
  occ <- tryCatch(
    CellMaNam::calc_occurrence(
      markers_data  = marker_ref,
      features_col  = "gene",
      cell_column   = "final_cluster"
    ),
    error = function(e)
      NULL
  )
  if (is.null(occ) ||
      !nrow(occ))
    return(setNames(character(0), character(0)))
  occ2 <- tryCatch(
    CellMaNam::select_top_occ(occ, top_n = top_n),
    error = function(e)
      NULL
  )
  if (is.null(occ2) ||
      !nrow(occ2))
    return(setNames(character(0), character(0)))
  ann_tbl <- tryCatch(
    CellMaNam::get_annotation(
      cell_markers = seurat_clusters_DEG %>% dplyr::select(cell_annotation = cluster, markers = gene),
      markers_occ  = occ2,
      max_genes    = nrow(
        Seurat::GetAssayData(seurat_obj, assay = cfg$base_assay, layer = "data")
      )
    ),
    error = function(e)
      NULL
  )
  if (is.null(ann_tbl) ||
      !nrow(ann_tbl))
    return(setNames(character(0), character(0)))
  cell_types <- tryCatch(
    CellMaNam::cell_typing(
      annotation_data = ann_tbl,
      hierarchy_data = NULL,
      p_val = p_val,
      level = level,
      hierarchy = FALSE
    ),
    error = function(e)
      NULL
  )
  if (is.null(cell_types) ||
      !nrow(cell_types))
    return(setNames(character(0), character(0)))
  need <- c("annotation", "full_names", "completed")
  if (!all(need %in% names(cell_types)))
    return(setNames(character(0), character(0)))
  df <- cell_types %>% dplyr::group_by(annotation) %>% dplyr::filter(completed == max(completed, na.rm =
                                                                                        TRUE)) %>% dplyr::ungroup() %>%
    dplyr::select(cluster = annotation, annotation = full_names)
  if (!nrow(df))
    return(setNames(character(0), character(0)))
  res <- df$annotation
  names(res) <- canon_cluster(df$cluster)
  res
}
cluster_expressed_bg <- function(seurat_obj, cluster_idents, min_cells = 5) {
  sct <- Seurat::GetAssayData(seurat_obj, assay = cfg$base_assay, layer = "data")
  g <- seurat_obj[[cluster_idents]][, 1] %>% canon_cluster()
  if (is.null(names(g)))
    names(g) <- colnames(seurat_obj)
  lapply(split(names(g), g), function(cells) {
    if (!length(cells))
      return(character(0))
    keep <- Matrix::rowSums(sct[, cells, drop = FALSE] > 0) >= min_cells
    rownames(sct)[keep]
  })
}
annotate_hypergeomX <- function(seurat_obj,
                                seurat_clusters_DEG,
                                marker_ref,
                                cluster_idents,
                                min_cells_bg = 5) {
  if (!nrow(seurat_clusters_DEG))
    return(setNames(character(0), character(0)))
  seurat_clusters_DEG$cluster <- canon_cluster(as.character(seurat_clusters_DEG$cluster))
  clusters <- unique(seurat_clusters_DEG$cluster)
  bg_by_cluster <- cluster_expressed_bg(seurat_obj, cluster_idents, min_cells =
                                          min_cells_bg)
  vals <- sapply(clusters, function(cl) {
    bg_genes <- unique(bg_by_cluster[[cl]])
    if (!length(bg_genes))
      return("Unknown")
    cl_genes <- intersect(unique(seurat_clusters_DEG$gene[seurat_clusters_DEG$cluster == cl]), bg_genes)
    if (!length(cl_genes))
      return("Unknown")
    df <- do.call(rbind, lapply(unique(marker_ref$final_cluster), function(ct) {
      ct_genes <- intersect(unique(marker_ref$gene[marker_ref$final_cluster == ct]), bg_genes)
      overlap <- length(intersect(cl_genes, ct_genes))
      m <- length(ct_genes)
      n <- length(bg_genes) - m
      k <- length(cl_genes)
      pval <- stats::phyper(overlap - 1, m, n, k, lower.tail = FALSE)
      data.frame(cell_type = ct, pval = pval)
    }))
    if (!nrow(df))
      return("Unknown")
    df$p_adj <- p.adjust(df$pval, method = "BH")
    as.character(df$cell_type[which.min(df$p_adj)])
  }, USE.NAMES = FALSE)
  names(vals) <- clusters
  vals
}
annotate_fgsea <- function(seurat_clusters_DEG,
                           marker_ref,
                           minSize = 3,
                           maxSize = 500) {
  if (!has_fgsea ||
      !nrow(seurat_clusters_DEG))
    return(setNames(character(0), character(0)))
  genesets <- split(marker_ref$gene, marker_ref$final_cluster)
  res <- lapply(split(seurat_clusters_DEG, seurat_clusters_DEG$cluster), function(df_cl) {
    lfc <- df_cl$avg_log2FC
    names(lfc) <- df_cl$gene
    if (!length(lfc))
      return("Unknown")
    lfc <- sort(tapply(lfc, names(lfc), max), decreasing = TRUE)
    gs <- lapply(genesets, function(v)
      intersect(v, names(lfc)))
    len <- sapply(gs, length)
    gs <- gs[len >= minSize & len <= maxSize]
    if (!length(gs))
      return("Unknown")
    gsr <- suppressWarnings(fgsea::fgsea(
      pathways = gs,
      stats = lfc,
      nperm = 2000
    ))
    if (!nrow(gsr))
      return("Unknown")
    as.character(gsr$pathway[order(gsr$padj, -abs(gsr$NES))][1])
  })
  labs <- unlist(res)
  names(labs) <- names(res)
  labs
}

# ---- UCell lean ----
fast_ucell_labels_lean <- function(seurat_obj,
                                   marker_ref,
                                   cluster_idents,
                                   min_genes   = cfg$ucell_min_genes,
                                   cells_per_cluster = cfg$ucell_cells_per_cluster,
                                   max_signatures    = cfg$ucell_max_signatures,
                                   ncores = cfg$ucell_ncores) {
  if (!has_ucell || !isTRUE(cfg$enable_ucell)) {
    message("[annot][ucell] UCell unavailable or disabled; skipping.")
    return(list(
      obj = seurat_obj,
      labels = setNames(character(0), character(0))
    ))
  }
  if (!cluster_idents %in% colnames(seurat_obj@meta.data)) {
    message("[annot][ucell] Cluster key '",
            cluster_idents,
            "' not found; skipping.")
    return(list(
      obj = seurat_obj,
      labels = setNames(character(0), character(0))
    ))
  }
  
  DefaultAssay(seurat_obj) <- cfg$base_assay
  bg <- rownames(Seurat::GetAssayData(seurat_obj, assay = cfg$base_assay, layer = "data"))
  if (!length(bg)) {
    message("[annot][ucell] No background genes in assay; skipping.")
    return(list(
      obj = seurat_obj,
      labels = setNames(character(0), character(0))
    ))
  }
  
  sigs <- split(marker_ref$gene, marker_ref$final_cluster)
  sigs <- lapply(sigs, function(v)
    intersect(unique(v), bg))
  sigs <- sigs[sapply(sigs, length) >= min_genes]
  if (!length(sigs)) {
    message(
      "[annot][ucell] No usable signatures after filtering (min_genes=",
      min_genes,
      "); skipping."
    )
    return(list(
      obj = seurat_obj,
      labels = setNames(character(0), character(0))
    ))
  }
  sig_len  <- sort(sapply(sigs, length), decreasing = TRUE)
  keep_sig <- names(sig_len)[seq_len(min(length(sig_len), max_signatures))]
  sigs     <- sigs[keep_sig]
  
  grp <- seurat_obj[[cluster_idents]][, 1] |> as.character() |> canon_cluster()
  if (length(grp) != ncol(seurat_obj)) {
    message("[annot][ucell] length(grp) != ncol(object); skipping.")
    return(list(
      obj = seurat_obj,
      labels = setNames(character(0), character(0))
    ))
  }
  cells_by <- split(colnames(seurat_obj), grp)
  subs <- unlist(lapply(cells_by, function(v)
    if (length(v) <= cells_per_cluster)
      v
    else
      sample(v, cells_per_cluster)), use.names = FALSE)
  if (!length(subs)) {
    message("[annot][ucell] Subsampling kept 0 cells; skipping.")
    return(list(
      obj = seurat_obj,
      labels = setNames(character(0), character(0))
    ))
  }
  
  sub <- subset(seurat_obj, cells = subs)
  pre_cols <- colnames(sub@meta.data)
  message(
    "[annot][ucell] Running UCell on ",
    length(subs),
    " cells across ",
    length(sigs),
    " signatures..."
  )
  sub <- UCell::AddModuleScore_UCell(sub,
                                     features = sigs,
                                     name = "U",
                                     ncores = ncores)
  
  post_cols <- setdiff(colnames(sub@meta.data), pre_cols)
  if (!length(post_cols))
    post_cols <- grep("^U[._-]", colnames(sub@meta.data), value = TRUE)
  if (!length(post_cols)) {
    message("[annot][ucell] No UCell score columns were created; skipping.")
    return(list(
      obj = seurat_obj,
      labels = setNames(character(0), character(0))
    ))
  }
  
  clean_uc <- function(x) {
    x <- gsub("^U[._-]*", "", x, perl = TRUE)
    x <- gsub("(?:[._-]*(?:UCell|U))?$", "", x, perl = TRUE)
    trim(x)
  }
  col2label <- setNames(clean_uc(post_cols), post_cols)
  
  sub_grp <- sub[[cluster_idents]][, 1] |> as.character() |> canon_cluster()
  df <- data.frame(cluster = sub_grp, sub@meta.data[, post_cols, drop = FALSE], check.names = FALSE)
  lab <- tapply(seq_len(nrow(df)), df$cluster, function(ix) {
    med <- suppressWarnings(apply(df[ix, post_cols, drop = FALSE], 2, stats::median, na.rm = TRUE))
    best_col <- names(which.max(med))
    col2label[[best_col]]
  })
  
  labs <- unlist(lab)
  labs <- setNames(as.character(labs), names(lab))
  list(obj = seurat_obj, labels = labs)
}

# Two-pass annotation orchestrator ---------------------------------------
annotate_all_methods <- function(obj,
                                 marker_ref,
                                 cluster_key_name = "cluster_key",
                                 prefix = "annot_detailed",
                                 limit_clusters = NULL,
                                 deg_prefix = "deg",
                                 degs_precomputed = NULL) {
  stopifnot(inherits(obj, "Seurat"))
  if (is.null(marker_ref) ||
      !nrow(marker_ref) ||
      !all(c("gene", "final_cluster") %in% names(marker_ref))) {
    stop("[annot] marker_ref must be a data.frame with columns: gene, final_cluster")
  }
  obj <- ensure_cluster_key(obj, cluster_key_name)
  
  old_id <- Idents(obj)
  on.exit(try(Idents(obj) <- old_id, silent = TRUE)
          , add = TRUE)
  Idents(obj) <- obj[[cluster_key_name]][, 1]
  clv <- Idents(obj)
  
  DefaultAssay(obj) <- cfg$base_assay
  X_full <- Seurat::GetAssayData(obj, assay = cfg$base_assay, layer = "data")
  if (!inherits(X_full, "dgCMatrix"))
    X_full <- methods::as(X_full, "dgCMatrix")
  bg <- rownames(X_full)
  
  curated_genes <- unique(as.character(marker_ref$gene))
  feats_cur <- intersect(curated_genes, bg)
  if (!length(feats_cur))
    stop("[annot] No curated markers found in assay background.")
  if (isTRUE(cfg$annot_features_max) &&
      length(feats_cur) > cfg$annot_features_max)
    feats_cur <- feats_cur[seq_len(cfg$annot_features_max)]
  
  obj_lean <- Seurat::DietSeurat(
    obj,
    assays = cfg$base_assay,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    dimreducs = character(),
    graphs = character(),
    features = feats_cur
  )
  DefaultAssay(obj_lean) <- cfg$base_assay
  Idents(obj_lean) <- obj_lean[[cluster_key_name]][, 1]
  
  if (!is.null(limit_clusters)) {
    limit_clusters <- canon_cluster(as.character(limit_clusters))
    grp_all <- as.character(Idents(obj_lean))
    keep_cells <- colnames(obj_lean)[grp_all %in% limit_clusters]
    if (!length(keep_cells)) {
      message("[annot] limit_clusters matched 0 cells; returning empty table.")
      tab_empty <- tibble::tibble(cluster = character(0))
      if (!prefix %in% colnames(obj@meta.data))
        obj[[prefix]] <- NA_character_
      if (!"final_consensus" %in% colnames(obj@meta.data))
        obj$final_consensus <- NA_character_
      try(ckpt_save("annot_table_last",
                    list(table = tab_empty, saved_at = Sys.time())), silent = TRUE)
      return(list(
        obj = obj,
        table = tab_empty,
        degs = tibble::tibble()
      ))
    }
    obj_lean <- subset(obj_lean, cells = keep_cells)
  }
  
  clusters <- Idents(obj_lean)
  lv <- levels(clusters)
  expected_clusters <- lv
  
  tab_placeholder <- tibble::tibble(cluster = canon_cluster(as.character(expected_clusters)))
  try(ckpt_save("annot_table_last",
                list(table = tab_placeholder, saved_at = Sys.time())), silent = TRUE)
  
  message("[annot] DEGs: preparing ...")
  degs <- degs_precomputed
  if (is.null(degs)) {
    if (isTRUE(cfg$annot_skip_deg)) {
      cached <- try(deg_ckpt_load(
        prefix = deg_prefix,
        obj = obj,
        group_col = cluster_key_name,
        features = feats_cur
      ),
      silent = TRUE)
      if (!inherits(cached, "try-error") &&
          !is.null(cached))
        degs <- cached
    } else {
      layer_or_slot <- pick_layer_arg()
      pieces <- lapply(expected_clusters, function(cl) {
        args <- list(
          object = obj_lean,
          ident.1 = cl,
          only.pos = TRUE,
          min.pct = 0.20,
          logfc.threshold = 0.25,
          test.use = "wilcox",
          verbose = FALSE,
          features = feats_cur,
          random.seed = cfg$seed
        )
        if (!is.null(cfg$deg_subsample_per_ident))
          args$max.cells.per.ident <- cfg$deg_subsample_per_ident
        if (layer_or_slot == "layer")
          args$layer <- "data"
        else
          args$slot <- "data"
        fm <- tryCatch(
          do.call(Seurat::FindMarkers, args),
          error = function(e)
            NULL
        )
        if (is.null(fm) || !nrow(fm))
          return(NULL)
        tibble::tibble(
          cluster    = canon_cluster(as.character(cl)),
          gene       = rownames(fm),
          avg_log2FC = suppressWarnings(
            if ("avg_log2FC" %in% names(fm))
              fm$avg_log2FC
            else if ("avg_logFC" %in% names(fm))
              fm$avg_logFC
            else if ("log2FC" %in% names(fm))
              fm$log2FC
            else
              NA_real_
          ),
          p_val_adj  = suppressWarnings(
            if ("p_val_adj"  %in% names(fm))
              fm$p_val_adj
            else if ("p_val.adj" %in% names(fm))
              fm$p_val.adj
            else if ("p_val" %in% names(fm))
              fm$p_val
            else
              1
          )
        )
      })
      degs <- dplyr::bind_rows(pieces)
      if (!is.null(degs) && nrow(degs)) {
        degs <- dplyr::filter(
          degs,
          is.finite(avg_log2FC),
          is.finite(p_val_adj),
          p_val_adj < 0.05,
          avg_log2FC > 0
        )
        try(deg_ckpt_save(
          prefix = deg_prefix,
          obj = obj,
          group_col = cluster_key_name,
          degs = degs,
          features = feats_cur
        ),
        silent = TRUE)
      }
    }
  }
  message(sprintf(
    "[annot] DEGs: %s",
    ifelse(
      is.null(degs) ||
        !nrow(degs),
      "none (empty table)",
      sprintf("n=%d", nrow(degs))
    )
  ))
  
  X <- Seurat::GetAssayData(obj_lean, assay = cfg$base_assay, layer = "data")
  stopifnot(inherits(X, "dgCMatrix"))
  avg_gc   <- avg_by_group_sparse(X, clusters)
  avg_gc_m <- as.matrix(avg_gc)
  
  obj_lean@misc$annot_feats   <- feats_cur
  obj_lean@misc$annot_avg_exp <- avg_gc_m
  obj_lean@misc$annot_degs    <- degs
  assign(".annot_obj_lean", obj_lean, envir = .GlobalEnv)
  
  out <- list()
  curated_genes <- unique(as.character(marker_ref$gene))
  out$avg_exp   <- .run_annot_method("avg_exp", function()
    annotate_avgexp_matrix(avg_gc_m, marker_ref, curated_genes = curated_genes), expect_clusters = expected_clusters)
  out$hypergeom <- .run_annot_method("hypergeom", function()
    annotate_hypergeom(degs %||% tibble::tibble(), marker_ref, rownames(X)), expect_clusters = expected_clusters)
  out$majority  <- .run_annot_method("majority", function()
    annotate_majority(degs %||% tibble::tibble(), marker_ref), expect_clusters = expected_clusters)
  out$logfc     <- .run_annot_method("logfc", function()
    annotate_logfc(degs %||% tibble::tibble(), marker_ref, use_padj_weight = TRUE), expect_clusters = expected_clusters)
  out$cellmanam <- .run_annot_method("cellmanam", function()
    annotate_cellmanam(obj, degs %||% tibble::tibble(), marker_ref), expect_clusters = expected_clusters)
  out$hypergeomX <- .run_annot_method("hypergeomX", function()
    annotate_hypergeomX(
      obj,
      degs %||% tibble::tibble(),
      marker_ref,
      cluster_idents = cluster_key_name,
      min_cells_bg = 5
    ), expect_clusters = expected_clusters)
  out$gsea      <- .run_annot_method("gsea", function()
    annotate_fgsea(degs %||% tibble::tibble(), marker_ref), expect_clusters = expected_clusters)
  
  uc <- fast_ucell_labels_lean(obj, marker_ref, cluster_idents = cluster_key_name)
  obj <- uc$obj
  out$ucell     <- .run_annot_method("ucell", function()
    uc$labels, expect_clusters = expected_clusters)
  
  build_tbl <- function(named_vec, nm) {
    if (!length(named_vec))
      return(NULL)
    tibble::tibble(cluster = names(named_vec), !!nm := unname(named_vec))
  }
  dfs <- purrr::compact(lapply(names(out), function(nm)
    build_tbl(out[[nm]], nm)))
  tab <- if (length(dfs))
    Reduce(function(x, y)
      dplyr::full_join(x, y, by = "cluster"), dfs)
  else
    tibble::tibble(cluster = character(0))
  if (!nrow(tab))
    tab <- tibble::tibble(cluster = canon_cluster(as.character(expected_clusters)))
  tab$cluster <- canon_cluster(tab$cluster)
  if (!is.null(limit_clusters))
    tab <- dplyr::filter(tab, cluster %in% expected_clusters)
  
  methods <- intersect(
    c(
      "avg_exp",
      "hypergeom",
      "majority",
      "logfc",
      "cellmanam",
      "hypergeomX",
      "gsea",
      "ucell"
    ),
    names(tab)
  )
  weights <- c(
    avg_exp = 1,
    hypergeom = 1,
    majority = 1,
    logfc = 1,
    cellmanam = 1,
    hypergeomX = 1,
    gsea = 1,
    ucell = 1
  )
  weights <- weights[methods]
  wvote <- function(row) {
    # methods + weights already defined above
    vals <- as.list(row[methods])
    labs <- vapply(vals, function(x) as.character(x[[1]]), "", USE.NAMES = FALSE)
    keep <- !is.na(labs) & labs != "Unknown"
    if (!any(keep)) return(NA_character_)
    # align weights to the methods that actually voted
    ms   <- methods[keep]
    labs <- labs[keep]
    # sum method weights per label
    sc <- tapply(weights[ms], labs, sum, simplify = TRUE)
    names(which.max(sc))
  }
  
  if (nrow(tab) &&
      length(methods))
    tab$final_consensus <- apply(tab[, methods, drop = FALSE], 1, wvote)
  else
    tab$final_consensus <- NA_character_
  
  if (nrow(tab) && length(methods)) {
    tab$agree <- apply(tab[, methods, drop = FALSE], 1, function(x) {
      v <- x[!is.na(x) & x != "Unknown"]
      if (!length(v))
        return(0)
      max(table(v)) / length(v)
    })
  } else
    tab$agree <- numeric(nrow(tab))
  
  grp <- obj[[cluster_key_name]][, 1] %>% as.character() %>% canon_cluster()
  names(grp) <- colnames(obj)
  init_meta <- function(n)
    rep(NA_character_, n)
  # existing block writes: annot_detailed_avg_exp, ..., annot_detailed_ucell
  for (m in methods) {
    colname <- paste0(prefix, "_", m)
    if (!colname %in% colnames(obj@meta.data))
      obj[[colname]] <- init_meta(ncol(obj))
    cmap <- setNames(tab[[m]], tab$cluster)
    obj[[colname]] <- apply_mapping(grp, cmap)
  }
  # also write UNPREFIXED aliases by request (result_obj$avg_exp, etc.)
  for (m in methods) {
    cmap <- setNames(tab[[m]], tab$cluster)
    obj[[m]] <- apply_mapping(grp, cmap)
  }
  
  if (!prefix %in% colnames(obj@meta.data))
    obj[[prefix]] <- init_meta(ncol(obj))
  if (nrow(tab))
    obj[[prefix]] <- apply_mapping(grp, setNames(tab$final_consensus, tab$cluster))
  if (!"final_consensus" %in% colnames(obj@meta.data))
    obj$final_consensus <- init_meta(ncol(obj))
  if (nrow(tab))
    obj$final_consensus <- apply_mapping(grp, setNames(tab$final_consensus, tab$cluster))
  
  tab <- dplyr::arrange(tab, cluster)
  try(ckpt_save("annot_table_last", list(table = tab, saved_at = Sys.time())), silent = TRUE)
  
  res <- list(
    obj = obj,
    table = tab,
    degs = degs %||% tibble::tibble()
  )
  class(res) <- c("annot_res", "list")
  return(res)
}

print_annotation_summary <- function(tab, n = NULL) {
  if (is.null(tab) ||
      !nrow(tab)) {
    message("[annot][summary] Empty annotation table.")
    return(invisible(NULL))
  }
  cols <- intersect(
    c(
      "cluster",
      "avg_exp",
      "hypergeom",
      "majority",
      "logfc",
      "cellmanam",
      "hypergeomX",
      "gsea",
      "ucell",
      "final_consensus",
      "agree"
    ),
    colnames(tab)
  )
  show <- tab[, cols, drop = FALSE] %>% dplyr::arrange(cluster)
  if (!is.null(n))
    show <- utils::head(show, n)
  print(show, row.names = FALSE)
  invisible(show)
}

write_results_xlsx <- function(wb, sheet_prefix, annot_tab, degs, anno_df) {
  addWorksheet(wb, paste0(sheet_prefix, "_annot"))
  writeData(wb, paste0(sheet_prefix, "_annot"), annot_tab)
  if (is.null(degs) ||
      !is.data.frame(degs) || !nrow(degs))
    return(invisible(NULL))
  if (!"cluster" %in% names(degs))
    stop("DEG table lacks 'cluster' column.")
  by_cluster <- split(degs, as.character(degs$cluster))
  adf_ok <- is.data.frame(anno_df) && nrow(anno_df) > 0
  if (adf_ok) {
    adf <- anno_df
    if (!"gene" %in% names(adf)) {
      alt <- intersect(names(adf),
                       c(
                         "Gene",
                         "GeneID",
                         "gene_id",
                         "Gene_name",
                         "GeneID_version"
                       ))
      if (length(alt))
        names(adf)[match(alt[1], names(adf))] <- "gene"
    }
    if (!"all_anno" %in% names(adf)) {
      alt <- intersect(
        c(
          "Annotation",
          "Uniprot_protein_name",
          "PFAM_domain_name",
          "NCBI_ID"
        ),
        names(adf)
      )
      adf$all_anno <- if (length(alt))
        do.call(paste, c(adf[alt], sep = "; "))
      else
        NA_character_
    }
    adf$gene <- trim(as.character(adf$gene))
    adf <- unique(adf[, c("gene", "all_anno")])
    adf$gene_nov <- sub("\\.[0-9]+$", "", adf$gene)
    adf_slim <- unique(adf[, c("gene", "gene_nov", "all_anno")])
  }
  
  norm_chr <- function(x)
    trim(as.character(x))
  for (nm in names(by_cluster)) {
    df <- by_cluster[[nm]]
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if (!"gene" %in% names(df))
      df$gene <- rownames(df)
    df$gene <- norm_chr(df$gene)
    if (adf_ok) {
      # 1) First join — be safe as well
      df1 <- dplyr::left_join(
        df,
        adf_slim[, c("gene","all_anno"), drop = FALSE],
        by = "gene"
      )
      df1$all_anno <- vapply(df1$all_anno, as_chr_collapse, "", USE.NAMES = FALSE)
      # 2) Only if any annotations are missing
      need <- is.na(df1$all_anno) | df1$all_anno == ""
      if (any(need)) {
        df1$gene_nov <- sub("\\.[0-9]+$", "", df1$gene)
        stopifnot(is.data.frame(df1[need, c("gene_nov"), drop = FALSE]))
        
        lk <- dplyr::left_join(
          df1[need, c("gene_nov"), drop = FALSE],                             # <- drop = FALSE
          unique(adf_slim[, c("gene_nov","all_anno"), drop = FALSE]),         # <- drop = FALSE
          by = "gene_nov"
        )
        
        df1$all_anno[need] <- lk$all_anno
        df1$gene_nov <- NULL
        df1$all_anno <- vapply(df1$all_anno, as_chr_collapse, "", USE.NAMES = FALSE)
      }
      df <- df1
      
      
    }
    sht <- sanitize_sheet(paste0(sheet_prefix, "_", nm))
    addWorksheet(wb, sht)
    writeData(wb, sht, df)
  }
}

# ---------------------------
# Cluster keys & subclustering
# ---------------------------

ensure_cluster_key <- function(obj, key = "cluster_key") {
  if (!(key %in% colnames(obj@meta.data))) {
    # build it from Idents if possible
    if (length(Idents(obj))) {
      obj[[key]] <- canon_cluster(as.character(Idents(obj)))
    } else {
      stop("'",
           key,
           "' not found in Seurat object meta.data and Idents() is empty.")
    }
  } else {
    obj[[key]] <- canon_cluster(as.character(obj[[key]][, 1]))
  }
  obj
}

.pretty_sizes <- function(tab_named)
  paste(sprintf("%d:%d", seq_along(tab_named), as.integer(tab_named)), collapse = ", ")
.min_allowed <- function(n, abs_min, prop_min, tiny_prop) {
  list(
    min_allowed   = max(abs_min, ceiling(prop_min  * n)),
    tiny_allowed  = max(abs_min, ceiling(tiny_prop * n)),
    abs_min       = abs_min,
    prop_min_calc = ceiling(prop_min  * n),
    tiny_min_calc = ceiling(tiny_prop * n)
  )
}
# --- replace the centroid merge helper with this hardened version ---
.merge_labels_by_centroid <- function(emb, lab, tiny_levels) {
  if (!length(tiny_levels))
    return(lab)
  lab <- as.character(lab)
  lab[is.na(lab)] <- "NA"                    # guard NAs
  tab <- table(lab)
  big_levels <- setdiff(names(tab), tiny_levels)
  if (!length(big_levels))
    return(lab)
  
  emb <- as.matrix(emb)
  stopifnot(nrow(emb) == length(lab))
  
  # centroids per level
  lab_f <- factor(lab)                # ensure factor
  cent  <- rowsum(emb, lab_f) / as.numeric(table(lab_f))
  cent <- as.matrix(cent)
  
  for (sm in tiny_levels) {
    # if tiny level vanished, skip
    if (!sm %in% rownames(cent) || !length(big_levels))
      next
    d <- sapply(big_levels, function(b) {
      v <- cent[sm, , drop = FALSE] - cent[b, , drop = FALSE]
      sum(v^2)
    })
    to <- names(which.min(d))
    lab[lab == sm] <- to
  }
  lab
}



# ---- CORE: robust subclustering for ONE parent  (FIXED) ---------------------
# ---- CORE: robust subclustering for ONE parent (never throws) ----
subcluster_one <- function(obj_in,
                           parent_key,
                           cluster_key_name = "cluster_key",
                           new_key_name     = "cluster_key_v2",
                           res_grid         = cfg$subclust_res,
                           npcs             = cfg$subclust_npcs,
                           k                = cfg$subclust_k,
                           min_abs          = cfg$min_child_abs,
                           min_prop         = cfg$min_child_prop,
                           tiny_prop        = cfg$subclust_tiny_prop,
                           collapse_tiny    = TRUE) {
  # --- helpers ---
  has_valid_graph <- function(o, gnm) {
    G <- o@graphs[[gnm]]
    if (is.null(G))
      return(FALSE)
    if (!methods::is(G, "dgCMatrix"))
      G <- methods::as(G, "dgCMatrix")
    Matrix::nnzero(G) > 0
  }
  pick_snn <- function(o) {
    cands <- names(o@graphs)
    if (!length(cands))
      return(NULL)
    snn <- cands[grepl("_snn", cands, fixed = TRUE)]
    if (!length(snn))
      snn <- cands[1]
    snn[1]
  }
  
  # --- find parent cells from PARENT labels only ---
  base_lab <- obj_in[[cluster_key_name]][, 1] |> as.character() |> canon_cluster()
  p_regex <- paste0("^", parent_key, "(\\.|$)")
  cells_parent <- colnames(obj_in)[grepl(p_regex, base_lab, perl = TRUE)]
  parent_n <- length(cells_parent)
  if (!parent_n) {
    message("[subclust] (", parent_key, ") no cells; skip.")
    return(obj_in)
  }
  
  # small parents: skip early
  if (parent_n < (cfg$sub_min_cells_for_split %||% 30L)) {
    message("[subclust] (keep ",
            parent_key,
            ") too few cells (n=",
            parent_n,
            "); skip.")
    return(obj_in)
  }
  
  thr <- .min_allowed(parent_n, min_abs, min_prop, tiny_prop)
  
  # subset
  DefaultAssay(obj_in) <- cfg$base_assay
  sub <- tryCatch(
    subset(obj_in, cells = cells_parent),
    error = function(e)
      NULL
  )
  if (is.null(sub) || ncol(sub) < 3L) {
    message("[subclust] (keep ", parent_key, ") subset invalid; skip.")
    return(obj_in)
  }
  
  # HVG + drop zero-variance
  if (length(VariableFeatures(sub)) < 200) {
    sub <- FindVariableFeatures(
      sub,
      selection.method = "vst",
      nfeatures = min(2000, nrow(sub)),
      verbose = FALSE
    )
  }
  feats <- intersect(VariableFeatures(sub), rownames(sub))
  if (!length(feats)) {
    message("[subclust] (keep ", parent_key, ") no HVGs; skip.")
    return(obj_in)
  }
  
  Xsub <- tryCatch(
    Seurat::GetAssayData(obj_in, assay = cfg$base_assay, layer = "data")[feats, cells_parent, drop = FALSE],
    error = function(e)
      NULL
  )
  if (is.null(Xsub)) {
    message("[subclust] (keep ", parent_key, ") no data layer; skip.")
    return(obj_in)
  }
  
  mu  <- Matrix::rowMeans(Xsub)
  mu2 <- Matrix::rowMeans(Xsub^2)
  vrs <- as.numeric(mu2 - mu^2)
  feats <- feats[vrs > 0]
  if (!length(feats)) {
    message("[subclust] (keep ", parent_key, ") no variable signal; skip.")
    return(obj_in)
  }
  
  # scaling
  sub <- tryCatch(
    ScaleData(sub, features = feats, verbose = FALSE),
    error = function(e)
      NULL
  )
  if (is.null(sub)) {
    message("[subclust] (keep ", parent_key, ") ScaleData failed; skip.")
    return(obj_in)
  }
  
  # PCA (exact SVD, bounded PCs)
  use_pcs_max <- max(2L, min(length(feats), parent_n - 2L))
  Xsc <- tryCatch(
    Seurat::GetAssayData(sub, assay = cfg$base_assay, layer = "scale.data")[feats, , drop = FALSE],
    error = function(e)
      NULL
  )
  if (is.null(Xsc) || nrow(Xsc) < 2L || ncol(Xsc) < 3L) {
    message("[subclust] (keep ",
            parent_key,
            ") scaled matrix too small; skip.")
    return(obj_in)
  }
  rank_est <- tryCatch(
    as.integer(Matrix::rankMatrix(t(Xsc))),
    error = function(e)
      NA_integer_
  )
  if (is.finite(rank_est))
    use_pcs_max <- max(2L, min(use_pcs_max, rank_est - 1L))
  use_pcs <- min(npcs, 50L, use_pcs_max)
  
  sub <- tryCatch(
    RunPCA(
      sub,
      features = rownames(Xsc),
      npcs = use_pcs,
      approx = FALSE,
      verbose = FALSE
    ),
    error = function(e)
      NULL
  )
  if (is.null(sub) || is.null(sub@reductions$pca)) {
    message("[subclust] (keep ", parent_key, ") PCA failed; skip.")
    return(obj_in)
  }
  
  emb <- tryCatch(
    Embeddings(sub, "pca"),
    error = function(e)
      NULL
  )
  if (is.null(emb) || ncol(emb) < 2L) {
    message("[subclust] (keep ",
            parent_key,
            ") PCA embeddings too small; skip.")
    return(obj_in)
  }
  emb <- as.matrix(emb)
  storage.mode(emb) <- "double"
  avail <- min(ncol(emb), use_pcs)
  
  # Neighbors with safe k
  k_eff <- max(5L, min(k, parent_n - 2L, max(10L, floor(parent_n * 0.25))))
  gname <- paste0(cfg$base_assay, "_snn_sub_k", k_eff)
  sub <- tryCatch(
    FindNeighbors(
      sub,
      dims = 1:avail,
      k.param = k_eff,
      graph.name = gname,
      verbose = FALSE
    ),
    error = function(e)
      NULL
  )
  if (is.null(sub)) {
    message("[subclust] (keep ",
            parent_key,
            ") FindNeighbors failed; skip.")
    return(obj_in)
  }
  
  # verify graph
  g_pick <- if (has_valid_graph(sub, gname))
    gname
  else
    pick_snn(sub)
  if (is.null(g_pick) || !has_valid_graph(sub, g_pick)) {
    message("[subclust] (keep ",
            parent_key,
            ") SNN graph empty/invalid; skip.")
    return(obj_in)
  }
  
  # scan resolutions, try multiple algorithms via safe_findclusters()
  found_multi <- FALSE
  last_ok <- NULL
  for (res in res_grid) {
    sub_try <- try(safe_findclusters(
      sub,
      graph.name = g_pick,
      resolution = res,
      seed = cfg$seed
    ),
    silent = TRUE)
    if (!inherits(sub_try, "try-error")) {
      cl <- as.character(Idents(sub_try))
      if (length(unique(cl)) > 1L) {
        sub <- sub_try
        found_multi <- TRUE
        break
      }
      last_ok <- sub_try
    }
  }
  # if never >1, keep latest valid object and skip
  if (!found_multi) {
    message("[subclust] (keep ",
            parent_key,
            ") only one child at all resolutions; skip.")
    return(obj_in)
  }
  
  # children sizes (pre-merge)
  cl <- as.character(Idents(sub))
  tab <- sort(table(cl), decreasing = TRUE)
  message(
    "[subclust] (",
    parent_key,
    ") parent_n=",
    parent_n,
    " | children: ",
    paste(sprintf("%d:%d", seq_along(tab), tab), collapse = ", "),
    " | min_allowed=",
    thr$min_allowed,
    " [abs=",
    thr$abs_min,
    ", prop=",
    min_prop,
    "→",
    thr$prop_min_calc,
    ", tiny=",
    tiny_prop,
    "→",
    thr$tiny_min_calc,
    "]"
  )
  
  # collapse tiny by centroid
  if (collapse_tiny && any(tab < thr$tiny_min_calc)) {
    tiny_lvls <- names(tab)[tab < thr$tiny_min_calc]
    cl <- .merge_labels_by_centroid(emb[, 1:avail, drop = FALSE], cl, tiny_lvls)
    Idents(sub) <- factor(cl)
    tab <- sort(table(cl), decreasing = TRUE)
    if (length(tab) <= 1L) {
      message("[subclust] (keep ",
              parent_key,
              ") tiny merge → single child; skip.")
      return(obj_in)
    }
  }
  
  # enforce min size
  if (any(tab < thr$min_allowed)) {
    small_lvls <- names(tab)[tab < thr$min_allowed]
    cl <- .merge_labels_by_centroid(emb[, 1:avail, drop = FALSE], cl, small_lvls)
    Idents(sub) <- factor(cl)
    tab <- sort(table(cl), decreasing = TRUE)
    if (length(tab) <= 1L) {
      message("[subclust] (keep ",
              parent_key,
              ") size merge → single child; skip.")
      return(obj_in)
    }
  }
  
  # final relabel gX.N
  ord <- names(tab)
  idx_map <- setNames(seq_along(ord), ord)
  child_nums   <- unname(idx_map[as.character(Idents(sub))])
  child_labels <- paste0(parent_key, ".", child_nums)
  names(child_labels) <- Cells(sub)
  
  # ensure v2 column exists as character
  if (!(new_key_name %in% colnames(obj_in@meta.data))) {
    obj_in@meta.data[[new_key_name]] <- as.character(obj_in[[cluster_key_name]][, 1])
  } else if (!is.character(obj_in@meta.data[[new_key_name]])) {
    obj_in@meta.data[[new_key_name]] <- as.character(obj_in@meta.data[[new_key_name]])
  }
  
  stopifnot(all(cells_parent %in% rownames(obj_in@meta.data)))
  # write values in the exact order of 'cells_parent'
  vals <- unname(child_labels[match(cells_parent, names(child_labels))])
  if (anyNA(vals)) {
    # if any NA slipped in, abort cleanly
    message("[subclust] (keep ",
            parent_key,
            ") child label mismatch; skip.")
    return(obj_in)
  }
  obj_in@meta.data[cells_parent, new_key_name] <- vals
  obj_in@meta.data[[new_key_name]] <- factor(obj_in@meta.data[[new_key_name]])
  
  message("[subclust] (",
          parent_key,
          ") accepted split → ",
          paste(sprintf("%d:%d", seq_along(tab), tab), collapse = ", "))
  obj_in
}




# --- STAGE: subclust driver -----------------------------------------------
auto_subcluster_suspicious <- function(obj_in,
                                       annot_table,
                                       cluster_key_name = "cluster_key",
                                       agree_cut = cfg$agree_cut) {
  stopifnot(inherits(obj_in, "Seurat"))
  obj_in <- ensure_cluster_key(obj_in, cluster_key_name)
  if (!all(c("cluster", "agree") %in% colnames(annot_table))) {
    message("[subclust] Annotation table lacks 'cluster' and/or 'agree'; skipping.")
    return(list(
      obj = obj_in,
      changed = FALSE,
      parents = character(0)
    ))
  }
  susp <- annot_table %>% dplyr::filter(agree <= agree_cut) %>% dplyr::arrange(agree) %>% dplyr::pull(cluster) %>% unique()
  nk_name <- paste0(cluster_key_name, "_v2")
  if (!length(susp)) {
    message("[subclust] No suspicious clusters (agree ≤ ",
            agree_cut,
            "); nothing to split.")
    return(list(
      obj = obj_in,
      changed = FALSE,
      parents = character(0)
    ))
  }
  message(
    "[subclust] Will subcluster ",
    length(susp),
    " parent(s): ",
    paste(utils::head(susp, 10), collapse = ", "),
    if (length(susp) > 10)
      " ..."
    else
      ""
  )
  obj_work <- obj_in
  accepted <- character(0)
  for (p in susp) {
    obj_try <- try(subcluster_one(
      obj_in = obj_work,
      parent_key = p,
      cluster_key_name = cluster_key_name,
      new_key_name = nk_name,
      collapse_tiny = TRUE
    ),
    silent = TRUE)
    if (inherits(obj_try, "try-error")) {
      msg <- tryCatch(
        conditionMessage(attr(obj_try, "condition")),
        error = function(e)
          "unknown error"
      )
      message("[subclust] (keep ", p, ") error: ", msg)
      next
    }
    if (nk_name %in% colnames(obj_try@meta.data)) {
      lab_now <- as.character(obj_try[[nk_name]][, 1])
      if (any(startsWith(lab_now, paste0(p, ".")), na.rm = TRUE)) {
        accepted <- c(accepted, p)
        obj_work <- obj_try
      } else {
        message("[subclust] (keep ", p, ") no new children written; skip.")
      }
    }
  }
  list(
    obj = obj_work,
    changed = length(accepted) > 0,
    parents = unique(accepted),
    new_key = nk_name
  )
}

needs_rerun_subclust <- function() {
  if (!ckpt_has("subclust"))
    return(TRUE)
  st <- ckpt_load("subclust")
  if (!isTRUE(st$changed))
    return(FALSE)
  nk <- st$new_key %||% "cluster_key_v2"
  obj2 <- st$obj
  is.null(obj2) || !(nk %in% colnames(obj2@meta.data))
}

# --- timing helpers ---
.stage_time_log_add <- function(stage,
                                t0,
                                t1,
                                ok = TRUE,
                                note = NULL) {
  rec <- data.frame(
    stage = stage,
    start = as.character(t0),
    end   = as.character(t1),
    elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")),
    ok = ok,
    note = if (is.null(note))
      ""
    else
      as.character(note),
    stringsAsFactors = FALSE
  )
  if (ckpt_has("timings")) {
    log <- ckpt_load("timings")
    log <- rbind(log, rec)
  } else
    log <- rec
  ckpt_save("timings", log)
  invisible(rec)
}
save_timings_to_xlsx <- function(wb, sheet_name = "timings") {
  if (!ckpt_has("timings"))
    return(invisible(NULL))
  log <- ckpt_load("timings")
  if (sheet_name %in% names(wb))
    removeWorksheet(wb, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, log)
  invisible(log)
}

# --- CONFIG KNOBS defaults (keep) ---
cfg$subclust_res        <- cfg$subclust_res        %||% c(0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0)
cfg$subclust_npcs       <- cfg$subclust_npcs       %||% 30L
cfg$subclust_k          <- cfg$subclust_k          %||% 30L
cfg$agree_cut           <- cfg$agree_cut           %||% 0.6
cfg$min_child_abs       <- cfg$min_child_abs       %||% 10L
cfg$min_child_prop      <- cfg$min_child_prop      %||% 0.005
cfg$subclust_tiny_prop  <- cfg$subclust_tiny_prop  %||% 0.015

# ---------------------------
# STAGES
# ---------------------------
stage_load <- function() {
  message("Stage load...")
  stopifnot(file.exists(cfg$paths$markers_xlsx))
  cell_markers <- openxlsx::read.xlsx(cfg$paths$markers_xlsx, sheet = 1)
  refs <- prepare_marker_ref(cell_markers, gene_col = "Markers_positive_SMESG")
  marker_ref_general  <- refs$general
  marker_ref_detailed <- refs$detailed
  marker_ref          <- marker_ref_detailed
  
  anno_df <- NULL
  if (!is.null(cfg$paths$anno_rdata) &&
      file.exists(cfg$paths$anno_rdata)) {
    load(cfg$paths$anno_rdata)
    nn <- ls()
    hit <- nn[sapply(nn, function(x)
      is.data.frame(get(x)) && "gene" %in% names(get(x)))]
    if (length(hit)) {
      anno_df <- get(hit[1])
      if (!"all_anno" %in% names(anno_df)) {
        alt <- intersect(
          c(
            "Annotation",
            "Uniprot_protein_name",
            "PFAM_domain_name",
            "NCBI_ID"
          ),
          names(anno_df)
        )
        anno_df$all_anno <- if (length(alt))
          do.call(paste, c(anno_df[alt], sep = "; "))
        else
          NA_character_
      }
      anno_df <- anno_df[, intersect(c("gene", "all_anno"), names(anno_df)), drop =
                           FALSE]
    }
  }
  if (!is.null(anno_df)) {
    # Ensure 'gene' and 'all_anno' exist and are plain character
    if (!"gene" %in% names(anno_df)) {
      alt <- intersect(names(anno_df),
                       c("Gene", "GeneID", "gene_id", "Gene_name", "GeneID_version"))
      if (length(alt)) names(anno_df)[match(alt[1], names(anno_df))] <- "gene"
    }
    if (!"all_anno" %in% names(anno_df)) {
      alt <- intersect(names(anno_df),
                       c("Annotation", "Uniprot_protein_name", "PFAM_domain_name", "NCBI_ID"))
      anno_df$all_anno <- if (length(alt)) do.call(paste, c(anno_df[alt], sep="; ")) else NA_character_
    }
    # <- NEW: flatten possible list column to character
    if (is.list(anno_df$all_anno)) {
      anno_df$all_anno <- vapply(anno_df$all_anno, as_chr_collapse, "", USE.NAMES = FALSE)
    } else {
      anno_df$all_anno <- as.character(anno_df$all_anno %||% "")
    }
    anno_df$gene <- trim(as.character(anno_df$gene))
    anno_df <- unique(anno_df[, c("gene", "all_anno")])
  }
  
  
  obj <- NULL
  if (!is.null(cfg$paths$matrix_rds) &&
      file.exists(cfg$paths$matrix_rds)) {
    message("Loading from matrix_rds...")
    mat <- readRDS(cfg$paths$matrix_rds)
    if (!inherits(mat, "dgCMatrix"))
      stop("matrix_rds must be dgCMatrix.")
    if (is.null(rownames(mat)))
      rownames(mat) <- paste0("g", seq_len(nrow(mat)))
    if (is.null(colnames(mat)))
      colnames(mat) <- paste0("c", seq_len(ncol(mat)))
    empty_counts <- new("dgCMatrix",
                        Dim = dim(mat),
                        Dimnames = dimnames(mat))
    obj <- CreateSeuratObject(
      counts = empty_counts,
      assay = cfg$base_assay,
      min.cells = 0,
      min.features = 0
    )
    if ("SetAssayData" %in% getNamespaceExports("SeuratObject")) {
      obj <- SeuratObject::SetAssayData(
        obj,
        assay = cfg$base_assay,
        layer = "data",
        new.data = mat
      )
      obj <- SeuratObject::SetAssayData(
        obj,
        assay = cfg$base_assay,
        layer = "counts",
        new.data = empty_counts
      )
    } else {
      obj[[cfg$base_assay]]@data   <- mat
      obj[[cfg$base_assay]]@counts <- empty_counts
    }
    DefaultAssay(obj) <- cfg$base_assay
  } else {
    stopifnot(file.exists(cfg$paths$seurat_rdata))
    message("Loading from seurat_rdata...")
    load(cfg$paths$seurat_rdata)
    if (exists("integrated_seurat_obj_annotated_new"))
      obj <- integrated_seurat_obj_annotated_new
    else if (exists("integrated_seurat_obj"))
      obj <- integrated_seurat_obj
    else {
      nn <- ls()
      cls <- sapply(nn, function(x)
        class(get(x, inherits = FALSE))[1])
      cand <- nn[grepl("Seurat", cls)]
      if (!length(cand))
        stop("No Seurat object found inside: ", cfg$paths$seurat_rdata)
      obj <- get(cand[1], inherits = FALSE)
    }
    DefaultAssay(obj) <- cfg$base_assay
  }
  
  obj_small <- diet_for_checkpoint(
    obj,
    keep_assay = cfg$base_assay,
    keep_reduc = character(),
    keep_graphs = character(),
    drop_counts = TRUE,
    drop_scale = TRUE
  )
  
  ckpt_save(
    "load",
    list(
      obj = obj_small,
      marker_ref_general  = marker_ref_general,
      marker_ref_detailed = marker_ref_detailed,
      marker_ref          = marker_ref,
      anno_df = anno_df
    )
  )
  invisible(TRUE)
}

ensure_scaled <- function(obj, assay, features, chunk = 4000L) {
  DefaultAssay(obj) <- assay
  present <- intersect(unique(features), rownames(Seurat::GetAssayData(
    obj, assay = assay, layer = "data"
  )))
  if (!length(present))
    return(obj)
  parts <- split(present, ceiling(seq_along(present) / max(1L, chunk)))
  for (fs in parts)
    obj <- ScaleData(
      object = obj,
      assay = assay,
      features = fs,
      do.center = TRUE,
      do.scale = TRUE,
      verbose = FALSE
    )
  obj
}

stage_pca <- function() {
  message("Stage PCA...")
  st  <- ckpt_load("load")
  obj <- st$obj
  DefaultAssay(obj) <- cfg$base_assay
  if (!length(VariableFeatures(obj)))
    obj <- FindVariableFeatures(
      object = obj,
      assay = cfg$base_assay,
      selection.method = "vst",
      nfeatures = 3000,
      verbose = FALSE
    )
  hvg <- VariableFeatures(obj)
  if (!length(hvg))
    stop("No variable features.")
  obj <- ensure_scaled(
    obj,
    assay = cfg$base_assay,
    features = hvg,
    chunk = 4000L
  )
  obj <- RunPCA(
    object = obj,
    features = hvg,
    npcs = cfg$max_pcs,
    reduction.name = cfg$pca_name,
    reduction.key = "PCu_",
    verbose = FALSE
  )
  stdev <- obj@reductions[[cfg$pca_name]]@stdev
  if (is.null(stdev) || !length(stdev))
    stop("PCA failed.")
  dims <- choose_pcs_by_knee(
    stdev,
    max_pcs = min(cfg$max_pcs, length(stdev)),
    variance_cut = cfg$variance_cut,
    smooth_k = cfg$knee_smooth,
    min_pcs = 25L
  )
  message("Chosen PCs: ", dims)
  ckpt_save(
    "pca",
    list(
      obj = obj,
      dims = dims,
      marker_ref_general  = st$marker_ref_general,
      marker_ref_detailed = st$marker_ref_detailed,
      marker_ref          = st$marker_ref,
      anno_df = st$anno_df
    )
  )
  invisible(TRUE)
}

stage_grid <- function() {
  message("Stage grid...")
  st  <- ckpt_load("pca")
  obj <- st$obj
  dims <- st$dims
  gs <- grid_search_clusters(
    obj,
    dims,
    target_n = cfg$target_n_clusters,
    k_grid   = cfg$k_grid %||% c(5L, 8L, 10L, 15L, 20L),
    res_init = cfg$res_init %||% 0.6,
    res_max = cfg$res_max %||% 10,
    max_steps = cfg$grid_max_steps %||% 20
  )
  obj_best   <- gs$obj
  best_stats <- gs$stats
  keep_graph <- gs$keep_graph
  message(
    sprintf(
      "Best: k=%s res=%s clusters=%s modularity=%s",
      best_stats$k,
      best_stats$resolution,
      best_stats$n_clusters,
      round(best_stats$modularity, 3)
    )
  )
  obj_small <- diet_for_checkpoint(obj_best,
                                   keep_reduc  = cfg$pca_name,
                                   keep_graphs = keep_graph)
  ckpt_save(
    "grid",
    list(
      obj        = obj_small,
      dims = dims,
      best_stats = best_stats,
      keep_graph = keep_graph,
      marker_ref_general  = st$marker_ref_general,
      marker_ref_detailed = st$marker_ref_detailed,
      marker_ref          = st$marker_ref,
      anno_df    = st$anno_df
    )
  )
  invisible(TRUE)
}

stage_umap <- function() {
  message("Stage UMAP...")
  st <- ckpt_load("grid")
  obj <- st$obj
  dims <- st$dims
  obj <- RunUMAP(
    obj,
    reduction = cfg$pca_name,
    dims = 1:dims,
    n.neighbors = 30,
    min.dist = 0.5,
    metric = "cosine",
    reduction.name = cfg$umap_name,
    seed.use = cfg$seed,
    verbose = FALSE
  )
  obj$cluster_key <- canon_cluster(as.character(Idents(obj)))
  obj_small <- diet_for_checkpoint(
    obj,
    keep_reduc = c(cfg$pca_name, cfg$umap_name),
    keep_graphs = st$keep_graph
  )
  ckpt_update("umap", st, obj = obj_small)
  invisible(TRUE)
}

stage_deg <- function() {
  message("Stage DEG...")
  st <- ckpt_load("umap")
  obj <- st$obj
  obj <- ensure_cluster_key(obj, "cluster_key")
  degs <- compute_degs_robust(obj, group_col = "cluster_key", features_whitelist = NULL)
  if (!nrow(degs))
    message("Stage DEG: no DEGs found.")
  ckpt_update("deg",
              st,
              obj = obj,
              degs = degs,
              group_col = "cluster_key")
  try(deg_ckpt_save(
    prefix = "deg",
    obj = obj,
    group_col = "cluster_key",
    degs = degs,
    features = NULL
  ),
  silent = TRUE)
  obj_small <- diet_for_checkpoint(
    obj,
    keep_reduc = c(cfg$pca_name, cfg$umap_name),
    keep_graphs = st$keep_graph
  )
  ckpt_update("umap", st, obj = obj_small)
  invisible(TRUE)
}

rehydrate_refs <- function(st = NULL) {
  try_ck <- function(tag)
    if (ckpt_has(tag))
      ckpt_load(tag)
  else
    NULL
  bins <- list(st,
               try_ck("umap"),
               try_ck("grid"),
               try_ck("pca"),
               try_ck("load"))
  bins <- bins[!vapply(bins, is.null, TRUE)]
  pick_df <- function(x)
    is.data.frame(x) &&
    nrow(x) > 0 && all(c("gene", "final_cluster") %in% names(x))
  pick_anno <- function(x)
    is.data.frame(x) && nrow(x) > 0 && "gene" %in% names(x)
  
  mr_gen <- mr_det <- mr <- anno_df <- NULL
  for (b in bins) {
    if (is.null(mr_gen) &&
        pick_df(b$marker_ref_general))
      mr_gen <- b$marker_ref_general
    if (is.null(mr_det) &&
        pick_df(b$marker_ref_detailed))
      mr_det <- b$marker_ref_detailed
    if (is.null(mr)     &&
        pick_df(b$marker_ref))
      mr     <- b$marker_ref
    if (is.null(anno_df) &&
        pick_anno(b$anno_df))
      anno_df <- b$anno_df
  }
  if (is.null(mr) || !pick_df(mr)) {
    stopifnot(file.exists(cfg$paths$markers_xlsx))
    cm  <- openxlsx::read.xlsx(cfg$paths$markers_xlsx, sheet = 1)
    refs <- prepare_marker_ref(cm, gene_col = "Markers_positive_SMESG")
    mr_gen <- refs$general
    mr_det <- refs$detailed
    mr <- mr_det
  }
  if (is.null(anno_df) &&
      !is.null(cfg$paths$anno_rdata) &&
      file.exists(cfg$paths$anno_rdata)) {
    load(cfg$paths$anno_rdata)
    nn <- ls()
    hit <- nn[sapply(nn, function(x)
      is.data.frame(get(x)) && "gene" %in% names(get(x)))]
    if (length(hit)) {
      anno_df <- get(hit[1])
      if (!"all_anno" %in% names(anno_df)) {
        alt <- intersect(
          c(
            "Annotation",
            "Uniprot_protein_name",
            "PFAM_domain_name",
            "NCBI_ID"
          ),
          names(anno_df)
        )
        anno_df$all_anno <- if (length(alt))
          do.call(paste, c(anno_df[alt], sep = "; "))
        else
          NA_character_
      }
      anno_df <- anno_df[, intersect(c("gene", "all_anno"), names(anno_df)), drop =
                           FALSE]
    }
  }
  list(
    marker_ref_general = mr_gen,
    marker_ref_detailed = mr_det,
    marker_ref = mr,
    anno_df = anno_df
  )
}

stage_annot1 <- function() {
  message("Stage annot1...")
  st <- ckpt_load("umap")
  refs <- rehydrate_refs(st)
  st$marker_ref_general  <- refs$marker_ref_general
  st$marker_ref_detailed <- refs$marker_ref_detailed
  st$marker_ref          <- refs$marker_ref
  st$anno_df             <- refs$anno_df
  ckpt_update(
    "umap",
    st,
    marker_ref_general = st$marker_ref_general,
    marker_ref_detailed = st$marker_ref_detailed,
    marker_ref = st$marker_ref,
    anno_df = st$anno_df
  )
  
  mr <- st$marker_ref %||% st$marker_ref_detailed %||% st$marker_ref_general
  obj <- st$obj
  
  degs_in <- NULL
  if (ckpt_has("deg")) {
    d <- ckpt_load("deg")
    if (is.list(d) && !is.null(d$degs))
      degs_in <- d$degs
  }
  if (is.null(degs_in)) {
    degs_in <- try(deg_ckpt_load(
      prefix = "deg",
      obj = obj,
      group_col = "cluster_key",
      features = NULL
    ),
    silent = TRUE)
    if (inherits(degs_in, "try-error"))
      degs_in <- NULL
  }
  mr <- st$marker_ref %||% st$marker_ref_detailed %||% st$marker_ref_general
  
  ann1 <- try(annotate_all_methods(
    obj,
    mr,
    cluster_key_name = "cluster_key",
    prefix = "annot_detailed",
    limit_clusters = NULL,
    deg_prefix = "deg",
    degs_precomputed = degs_in
  ),
  silent = TRUE)
  
  good <- is.list(ann1) &&
    !is.null(ann1$obj) && !is.null(ann1$table)
  if (!good) {
    msg <- if (inherits(ann1, "try-error"))
      as.character(attr(ann1, "condition"))
    else
      "non-list return"
    message("[annot1] annotate_all_methods failed: ", msg)
    message("[annot1] Falling back to a minimal round1 table so the pipeline can proceed.")
    obj <- ensure_cluster_key(obj, "cluster_key")
    parents <- sort(unique(canon_cluster(as.character(
      obj$cluster_key
    ))))
    tab1 <- tibble::tibble(cluster = parents,
                           final_consensus = NA_character_,
                           agree = 0)
    try(ckpt_save("annot_table_last", list(table = tab1, saved_at = Sys.time())), silent = TRUE)
    wb <- wb_load_or_new(cfg$paths$out_xlsx)
    write_results_xlsx(
      wb,
      sheet_prefix = "round1",
      annot_tab = tab1,
      degs = if (isTRUE(cfg$write_round1_degs) &&
                 !is.null(degs_in) && nrow(degs_in))
        degs_in
      else
        NULL,
      anno_df = st$anno_df
    )
    save_timings_to_xlsx(wb)
    wb_save(wb, cfg$paths$out_xlsx)
    obj_small <- diet_for_checkpoint(
      obj,
      keep_reduc  = c(cfg$pca_name, cfg$umap_name),
      keep_graphs = st$keep_graph
    )
    ckpt_update("annot1", st, obj = obj_small, tab1 = tab1)
    return(invisible(TRUE))
  }
  
  obj  <- ann1$obj
  tab1 <- ann1$table %>% dplyr::arrange(cluster)
  message("[annot][summary] Round1 (detailed) first 20 rows:")
  invisible(print_annotation_summary(tab1, n = 20))
  try(ckpt_save("annot1_table", list(tab = tab1, saved_at = Sys.time())), silent = TRUE)
  
  wb <- wb_load_or_new(cfg$paths$out_xlsx)
  write_results_xlsx(
    wb,
    sheet_prefix = "round1",
    annot_tab = tab1,
    degs = if (isTRUE(cfg$write_round1_degs))
      ann1$degs
    else
      NULL,
    anno_df = st$anno_df
  )
  save_timings_to_xlsx(wb)
  wb_save(wb, cfg$paths$out_xlsx)
  
  obj_small <- diet_for_checkpoint(
    obj,
    keep_reduc  = c(cfg$pca_name, cfg$umap_name),
    keep_graphs = st$keep_graph
  )
  ckpt_update("annot1", st, obj = obj_small, tab1 = tab1)
  invisible(TRUE)
}

auto_split_large <- function(obj,
                             key        = "cluster_key_v2",
                             max_cells  = 1200L,
                             max_passes = 2L) {
  for (pass in seq_len(max_passes)) {
    lab <- as.character(obj[[key]][,1])
    tt  <- sort(table(lab), decreasing = TRUE)
    big <- names(tt)[tt > max_cells]
    if (!length(big)) {
      message("[autosplit] pass ", pass, ": nothing above ", max_cells, " cells.")
      break
    }
    message("[autosplit] pass ", pass, ": splitting ", length(big), " parents: ",
            paste(head(big, 10), collapse=", "), if (length(big)>10) " ..." else "")
    for (p in big) {
      obj <- subcluster_one(
        obj_in = obj,
        parent_key = p,
        cluster_key_name = key,
        new_key_name     = key,   # in-place refinement of the same key
        collapse_tiny    = TRUE
      )
    }
  }
  obj$cluster_key_final <- obj[[key]][,1]
  obj
}

stage_subclust <- function() {
  message("Stage subcluster...")
  st <- ckpt_load("annot1")
  if (is.null(st))
    stop("[subclust] Missing 'annot1' checkpoint. Run annot1 first.")
  obj  <- st$obj
  tab1 <- st$tab1
  obj <- ensure_cluster_key(obj, "cluster_key")
  Idents(obj) <- obj$cluster_key
  # if v2 exists but has no children (no dots), drop it to avoid confusion
  if ("cluster_key_v2" %in% colnames(obj@meta.data) &&
      !any(grepl("\\.", obj$cluster_key_v2))) {
    obj$cluster_key_v2 <- NULL
  }
  
  
  susp <- tryCatch(
    tab1 %>% dplyr::filter(!is.na(agree) &
                             agree <= cfg$agree_cut) %>% dplyr::pull(cluster) %>% unique() %>% as.character(),
    error = function(e)
      character(0)
  )
  if (!length(susp)) {
    message("[subclust] No parents below agree_cut (",
            cfg$agree_cut,
            "). Nothing to split.")
    ckpt_save("subclust", c(st, list(obj = obj, changed = FALSE)))
    return(invisible(TRUE))
  }
  message(
    "[subclust] Will subcluster ",
    length(susp),
    " parent(s): ",
    paste(head(susp, 10), collapse = ", "),
    if (length(susp) > 10)
      " ..."
    else
      ""
  )
  
  res <- auto_subcluster_suspicious(
    obj,
    annot_table = tab1,
    cluster_key_name = "cluster_key",
    agree_cut = cfg$agree_cut
  )
  obj <- res$obj
  
  if (!isTRUE(res$changed)) {
    message("[subclust] No accepted splits (",
            length(res$parents),
            " accepted).")
    ckpt_save("subclust", c(st, list(obj = obj, changed = FALSE)))
    return(invisible(TRUE))
  }
  
  gd <- tryCatch(
    ckpt_load("grid"),
    error = function(e)
      NULL
  )
  dims_umap <- if (!is.null(gd) && length(gd$dims))
    gd$dims
  else
    30L
  obj <- RunUMAP(
    obj,
    reduction = cfg$pca_name,
    dims = 1:dims_umap,
    n.neighbors = 30,
    min.dist = 0.5,
    metric = "cosine",
    reduction.name = paste0(cfg$umap_name, ".v2"),
    seed.use = cfg$seed,
    verbose = FALSE
  )
  
  #ckpt_save("subclust", c(st, list(obj = obj, parents = unique(res$parents), new_key = res$new_key, changed = TRUE)))
  st$obj     <- obj
  st$parents <- unique(res$parents)
  st$new_key <- res$new_key
  st$changed <- TRUE
  tot <- ncol(st$obj)
  max_cells_auto <- max(800L, round(0.0125 * tot))  #  ~1.25% of total cells, but at least 800
  st$obj <- auto_split_large(st$obj, key = "cluster_key_v2", max_cells = max_cells_auto, max_passes = 2L)
  ckpt_save("subclust", st)
  invisible(TRUE)
}




stage_annot2_and_final <- function() {
  message("Stage annot2/final...")
  st1 <- ckpt_load("annot1")
  stopifnot(!is.null(st1))
  obj  <- st1$obj
  tab1 <- st1$tab1
  
  # Use subclustered key if available/meaningful
  sc <- tryCatch(
    ckpt_load("subclust"),
    error = function(e)
      NULL
  )
  use_v2 <- FALSE
  if (!is.null(sc) && !is.null(sc$obj)) {
    if ("cluster_key_v2" %in% colnames(sc$obj@meta.data)) {
      v2 <- as.character(sc$obj$cluster_key_v2)
      use_v2 <- any(grepl("\\.", v2), na.rm = TRUE) ||
        (length(unique(na.omit(v2))) > length(unique(as.character(
          sc$obj$cluster_key
        ))))
      if (use_v2)
        obj <- sc$obj
    }
  }
  final_key <- if (use_v2)
    "cluster_key_v2"
  else
    "cluster_key"
  message(
    if (use_v2)
      "[annot2] Using subclustered key (cluster_key_v2)."
    else
      "[annot2] Using round1 key (cluster_key)."
  )
  obj$cluster_key_final <- obj[[final_key]][, 1]
  
  # Rehydrate marker refs / annotations
  refs <- rehydrate_refs()
  mr   <- refs$marker_ref %||% refs$marker_ref_detailed %||% refs$marker_ref_general
  
  # NEW: compute DEGs for the FINAL grouping and cache them (prefix 'deg_final')
  degs_final <- compute_degs_robust(obj, group_col = "cluster_key_final", features_whitelist = NULL)
  if (nrow(degs_final)) {
    try(deg_ckpt_save(
      prefix = "deg_final",
      obj = obj,
      group_col = "cluster_key_final",
      degs = degs_final,
      features = NULL
    ),
    silent = TRUE)
  } else {
    message(
      "[annot2] WARNING: no DEGs passed filters for final key; DEG-based methods may return 0 labels."
    )
  }
  
  # Re-annotate using final key; feed the DEGs we just computed
  ann2 <- annotate_all_methods(
    obj,
    mr,
    cluster_key_name = "cluster_key_final",
    prefix           = "annot_final",
    limit_clusters   = NULL,
    deg_prefix       = "deg_final",
    degs_precomputed = if (nrow(degs_final))
      degs_final
    else
      NULL
  )
  obj  <- ann2$obj
  tab2 <- ann2$table %>% dplyr::arrange(cluster)
  
  # Attach parent info from round-1 so it’s visible in Excel
  # Force data.frame/tibble types before joining (prevents class='character' surprises)
  # ---- SAFER: attach parent info without dplyr joins ----
  # Ensure tab2 is a data.frame with a 'cluster' column
  tab2 <- ann2$table
  if (!is.data.frame(tab2))
    tab2 <- as.data.frame(tab2, stringsAsFactors = FALSE)
  if (!"cluster" %in% names(tab2)) {
    if (!is.null(rownames(tab2))) {
      tab2$cluster <- rownames(tab2)
    } else {
      stop("[annot2] 'tab2' lacks a 'cluster' column and rownames; cannot augment.")
    }
  }
  tab2$cluster <- as.character(tab2$cluster)
  tab2 <- tab2[order(tab2$cluster), , drop = FALSE]
  
  # Map each final label to its parent
  df_labels <- data.frame(cluster = sort(unique(as.character(
    obj$cluster_key_final
  ))), stringsAsFactors = FALSE)
  df_labels$parent <- sub("\\..*$", "", df_labels$cluster)
  
  # Slim round1 table (parent consensus/agree)
  tab1_slim <- as.data.frame(st1$tab1, stringsAsFactors = FALSE)
  tab1_slim <- tab1_slim[, c("cluster", "final_consensus", "agree")]
  names(tab1_slim) <- c("parent", "parent_consensus", "parent_agree")
  tab1_slim$parent <- as.character(tab1_slim$parent)
  message("tab2 class: ", paste(class(ann2$table), collapse = ", "))
  message("tab1 class: ", paste(class(st1$tab1), collapse = ", "))
  
  # Base merges (keeps left order; no S3 generics involved)
  tab2_tmp <- merge(tab2,
                    df_labels,
                    by = "cluster",
                    all.x = TRUE,
                    sort = FALSE)
  tab2_aug <- merge(tab2_tmp,
                    tab1_slim,
                    by = "parent",
                    all.x = TRUE,
                    sort = FALSE)
  # --------------------------------------------------------
  
  # Write Excel: final_annot (augmented) + per-cluster DEG sheets
  wb <- wb_load_or_new(cfg$paths$out_xlsx)
  write_results_xlsx(
    wb,
    sheet_prefix = "final",
    annot_tab    = tab2_aug,
    degs         = ann2$degs,
    # will create final_<cluster> sheets now that DEGs exist
    anno_df      = refs$anno_df
  )
  save_timings_to_xlsx(wb)
  wb_save(wb, cfg$paths$out_xlsx)
  
  # Save checkpoint
  obj_small <- diet_for_checkpoint(obj, keep_reduc = c(cfg$pca_name, cfg$umap_name, paste0(cfg$umap_name, ".v2")))
  ckpt_save(
    "annot2_and_final",
    list(
      obj          = obj_small,
      final_key    = "cluster_key_final",
      final_labels = tab2_aug
    )
  )
  invisible(TRUE)
}


# Orchestrator ------------------------------------------------------------
auto_run <- function(resume = TRUE,
                     stop_after = NULL) {
  stages <- c("load",
              "pca",
              "grid",
              "umap",
              "deg",
              "annot1",
              "subclust",
              "annot2_and_final")
  fns <- list(
    load = stage_load,
    pca = stage_pca,
    grid = stage_grid,
    umap = stage_umap,
    deg = stage_deg,
    annot1 = stage_annot1,
    subclust = stage_subclust,
    annot2_and_final = stage_annot2_and_final
  )
  for (s in stages) {
    skip <- resume && ckpt_has(s)
    if (s == "subclust" &&
        resume && ckpt_has("subclust") && needs_rerun_subclust()) {
      message("[auto_run] subclust checkpoint is stale; re-running subclust.")
      skip <- FALSE
    }
    if (skip) {
      message("Skipping stage ", s, " (found checkpoint).")
      next
    }
    message(">>> Running stage ", s)
    t0 <- Sys.time()
    ok <- TRUE
    note <- NULL
    tryCatch({
      fns[[s]]()
    }, error = function(e) {
      ok <- FALSE
      note <- conditionMessage(e)
      stop(e)
    }, finally = {
      t1 <- Sys.time()
      .stage_time_log_add(
        stage = s,
        t0 = t0,
        t1 = t1,
        ok = ok,
        note = note
      )
      message(sprintf("[timing] %s: %.1f sec", s, as.numeric(difftime(t1, t0, units =
                                                                        "secs"))))
    })
    if (!is.null(stop_after) &&
        identical(s, stop_after)) {
      message("Stop requested after stage: ", s)
      break
    }
  }
  invisible(TRUE)
}

# ---------------------------
# RUN
# ---------------------------
CKPT_DIR  <- cfg$ckpt_dir
ckpt_file <- function(stage)
  file.path(CKPT_DIR, paste0("auto_annot_ckpt_", stage, ".rds"))

setwd("G:/PhD_final/sncRNA")

# Clear only late-stage checkpoints if needed
# unlink(ckpt_path_qs("annot1"),          force = TRUE); unlink(ckpt_path_rds("annot1"),          force = TRUE)
# unlink(ckpt_path_qs("subclust"),        force = TRUE); unlink(ckpt_path_rds("subclust"),        force = TRUE)
# unlink(ckpt_path_qs("annot2_and_final"), force = TRUE)
# unlink(ckpt_path_rds("annot2_and_final"), force = TRUE)
# unlink(ckpt_path_qs("final"), force = TRUE)
# unlink(ckpt_path_rds("final"), force = TRUE)

auto_run(resume = TRUE)

# Convenience diagnostics:
# diag <- if (ckpt_has("grid"))
#   ckpt_load("grid")$obj@misc$grid_diag
# else
#   NULL

# Convenience: pull final object (saved in annot2_and_final)
if (ckpt_has("annot2_and_final")) {
  result_obj <- ckpt_load("annot2_and_final")$obj
  message("Final object available as `result_obj` (labels in '$cluster_key_final').")
}
if (exists(".annot_obj_lean", envir = .GlobalEnv)) {
  obj_lean <- get(".annot_obj_lean", envir = .GlobalEnv)
  X <- Seurat::GetAssayData(obj_lean, assay = cfg$base_assay, layer = "data")
  print(class(X))
  print(dim(X))
  print(object.size(X))
  str(obj_lean@misc$annot_avg_exp)
}
st <- ckpt_load("annot1")
v <- as.character(st$obj$cluster_key)   # st <- ckpt_load("annot1")
names(v)
sc <- ckpt_load("subclust")
table(names(sc))["obj"]  # will be 2

#result_obj
#result_obj$final_consensus
# 1) Confirm final key has children
table(grepl("\\.", result_obj$cluster_key_final))
head(sort(unique(
  as.character(result_obj$cluster_key_final)
)))

# 2) Confirm per-method unprefixed columns exist
head(colnames(result_obj@meta.data)[colnames(result_obj@meta.data) %in%
                                      c(
                                        "avg_exp",
                                        "hypergeom",
                                        "majority",
                                        "logfc",
                                        "cellmanam",
                                        "hypergeomX",
                                        "gsea",
                                        "ucell",
                                        "final_consensus",
                                        "final_agree"
                                      )])

# 3) Excel: look for a 'final_annot' sheet and per-cluster 'final_<cluster>' sheets (with DEGs)
table(grepl("\\.", result_obj$cluster_key_final))           # subclusters present
head(result_obj@meta.data$avg_exp)                           # unprefixed methods exist
"final_annot" %in% openxlsx::getSheetNames(cfg$paths$out_xlsx)
st <- ckpt_load("annot1")
obj <- st$obj
Idents(obj) <- obj$cluster_key
obj2 <- subcluster_one(
  obj_in = obj,
  parent_key = "g60",
  cluster_key_name = "cluster_key",
  new_key_name = "cluster_key_v2"
)
# Should not error; either split or skip cleanly:
table(startsWith(as.character(obj2$cluster_key_v2), "g60."))
unique(result_obj$annot_final)
#result_obj$cluster_key_final
warnings()
head(result_obj$cluster_key_final)
unique(result_obj$cluster_key_final)
unique(result_obj$cluster_key_v2)
unique(result_obj$final_consensus)
