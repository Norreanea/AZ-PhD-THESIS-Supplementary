###############################################################################
# small RNA activity / enrichment in scRNA-seq (no Wilcoxon; no UCell)
#
# Core idea
# 1) Build sRNA→target gene sets by seed scanning against UTR (± CDS).
# 2) Score each cell by a control-matched module score for the target set.
#    - This reduces “cell type baseline” bias that often makes one population dominate.
# 3) Summarize per (celltype × timepoint × genotype) and test genotype effects
#    using a permutation test (directional, based on bulk sRNA change).
#
# Notes / constraints
# - You have 1 library per (condition × timepoint), so genotype p-values are
#   exploratory (cells are not true biological replicates). Treat as prioritization.
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(1)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(Biostrings)
  library(data.table)
  library(stringr)
  library("xlsx")
  library(ggplot2)
})
# ---------------------- sRNA sequence extraction -----------------------------
extract_seq_from_id <- function(id) {
  # pull a plausible nucleotide string; supports N and U
  s <- stringr::str_extract(id, "[ACGTUNacgtun]{15,100}")
  if (is.na(s)) return(NA_character_)
  s <- toupper(s)
  s <- chartr("U", "T", s)
  s
}
# ------------------------------- load data -----------------------------------
load("G:/PhD_final/result_obj_new.RData")                  # result_obj: Seurat obj
#load("G:/PhD_final/tRNA_miRNA_selected_raw_counts.RData")  # tRNA_miRNA_selected_raw_counts
load("G:/PhD_final/all_stringtie_selected.RData")          # all_stringtie_selected
load("D:/scRNA-seq/tRF_motif/cds_tx_seq.RData")            # cds_tx_seq
load("D:/scRNA-seq/tRF_motif/tx2gene.RData")               # tx2gene
load("D:/scRNA-seq/tRF_motif/utr_tx_seq.RData")            # utr_tx_seq
load("D:/scRNA-seq/AZ_final_obj/filtered_DEG_abr_new.RData")# filtered_DE (optional)
load("G:/PhD_final/final_bulk_DGE.RData")                  # final_test_DGE (optional)
#load("G:/PhD_final/tables/bulk_dir_tbl_LFC057.RData")      # bulk_dir_tbl (sRNA DE)

#all RNA type accumulation
not_rRNA=read.xlsx("D:/Elac2/final_results/tables/DGE_other_than_rRF_snRNA_filtered_new.xlsx",sheetIndex=1)
unique(not_rRNA$RNA_type)
bulk_dir_tbl <- not_rRNA
bulk_dir_tbl <- bulk_dir_tbl[bulk_dir_tbl$set=="Elac_vs_WT_dpa3" & bulk_dir_tbl$RNA_type %in% c("miRNAs","piRNAs","tRFs"),]
bulk_dir_tbl <- bulk_dir_tbl[,c("snRNA_type","Sequence","log2FoldChange")]
colnames(bulk_dir_tbl) <- c("snRNA_type","Sequence","bulk_log2FC_3dpa")
bulk_dir_tbl$Sequence <- chartr("U", "T", bulk_dir_tbl$Sequence)
bulk_dir_tbl$sRNA <- paste(bulk_dir_tbl$Sequence,bulk_dir_tbl$snRNA_type,sep=" ")

# ----------------------------- user parameters --------------------------------
OUT_DIR <- "G:/PhD_final/tables/srna_activity_sc"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Seurat parsing
ASSAY_USE <- "RNA"              # or "SCT" if you want SCT residual-like data
SLOT_USE  <- "data"             # "data" = log-normalized; "counts" only if you know what you're doing
MIN_CELLS_PER_GROUP <- 20       # per (celltype,timepoint,genotype)

# Which genotypes to compare
GENO_KEEP <- c("WT", "ELAC", "GFP")    # you can include "GFP" if needed

# Timepoints to analyze (NULL = all detected)
#TIMEPOINTS_TO_USE <- c(0, 16, 24, 72)
TIMEPOINTS_TO_USE <- 72
# Target scanning
USE_CDS_IN_SCAN <- TRUE
CDS_WEIGHT_MULT <- 0.5          # downweight CDS sites vs UTR
MIN_GENE_SCORE  <- 3            # minimum weighted site score to keep a gene
TOP_N_GENES     <- 500          # cap target set size (keeps scoring stable)

# Module score (control-matched)
N_BINS_EXPR     <- 24           # expression bins for control matching
CTRL_PER_TARGET <- 20           # controls sampled per target gene
MIN_TARGETS_FOR_SCORE <- 15

# Permutation test (run only on top strata per sRNA to keep runtime sane)
DO_PERMUTATION  <- TRUE
N_PERM          <- 2000
TEST_TOP_STRATA <- 15           # per (sRNA,model): compute perm p only for top |delta| strata

# If bulk_dir_tbl has expected_target_change (DOWN_in_ELAC/UP_in_ELAC) use it; else derive from bulk_log2FC_3dpa
# expected_target_change = "DOWN_in_ELAC" means sRNA is UP in ELAC (targets expected DOWN) => activity expected HIGHER in ELAC.

TARGET_MODELS <- c("miRNA_canonical", "piRNA_extended", "off1_7mer", "off2_7mer", "off3_7mer")

# Cache (target scanning is expensive)
TARGET_CACHE_RDS <- file.path(OUT_DIR, "targets_cache_corr.rds")

# ----------------------------- sanity checks ---------------------------------
stopifnot(exists("result_obj"))
seu <- result_obj
Seurat::DefaultAssay(seu) <- ASSAY_USE

stopifnot(exists("bulk_dir_tbl"))
bulk_dir_tbl <- as.data.table(bulk_dir_tbl)
stopifnot("sRNA" %in% names(bulk_dir_tbl))

if (!inherits(utr_tx_seq, "DNAStringSet")) utr_tx_seq <- Biostrings::DNAStringSet(utr_tx_seq)
if (!inherits(cds_tx_seq, "DNAStringSet")) cds_tx_seq <- Biostrings::DNAStringSet(cds_tx_seq)

# tx2gene can be a named vector or a 2-col df
as_tx2gene_named <- function(tx2gene_obj) {
  if (is.vector(tx2gene_obj) && !is.null(names(tx2gene_obj))) return(tx2gene_obj)
  if (is.data.frame(tx2gene_obj)) {
    cn <- tolower(colnames(tx2gene_obj))
    tx_col   <- colnames(tx2gene_obj)[match(TRUE, cn %in% c("tx","transcript","transcript_id","tx_id"))]
    gene_col <- colnames(tx2gene_obj)[match(TRUE, cn %in% c("gene","gene_id","geneid"))]
    if (is.na(tx_col) || is.na(gene_col)) stop("tx2gene needs transcript and gene columns.")
    v <- as.character(tx2gene_obj[[gene_col]])
    names(v) <- as.character(tx2gene_obj[[tx_col]])
    return(v)
  }
  stop("Unrecognized tx2gene format.")
}
tx2gene_map <- as_tx2gene_named(tx2gene)

# ------------------------- metadata column discovery --------------------------
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[[1]])
  NA_character_
}

meta <- seu@meta.data

CELLTYPE_COL <- pick_col(meta, c("final_population","celltype_use","celltype","CellType","seurat_clusters"))
COND_COL     <- pick_col(meta, c("condition_correct","sc_condition_full","condition","orig.ident"))

unique(meta$final_population)
if (is.na(CELLTYPE_COL) || is.na(COND_COL)) {
  stop("Could not find metadata columns for celltype and condition. Update CELLTYPE_COL / COND_COL candidates.")
}

# Parse strings like "WT13S", "ELAC24S", "GFP0S" -> genotype + timepoint
parse_condition <- function(x) {
  x <- as.character(x)
  base <- stringr::str_extract(x, "^[A-Za-z]+")
  tp   <- suppressWarnings(as.integer(stringr::str_extract(x, "[0-9]+")))
  data.table(genotype = base, timepoint = tp, condition_full = x)
}

parsed <- parse_condition(meta[[COND_COL]])
seu$genotype <- parsed$genotype
seu$timepoint <- parsed$timepoint
seu$celltype_use <- as.character(meta[[CELLTYPE_COL]])

# keep timepoints
tp_present <- sort(unique(seu$timepoint[!is.na(seu$timepoint)]))
TP_USE <- if (is.null(TIMEPOINTS_TO_USE)) tp_present else intersect(tp_present, TIMEPOINTS_TO_USE)
if (length(TP_USE) == 0) stop("No requested timepoints found in Seurat metadata.")

# subset
cells_keep <- which(seu$genotype %in% GENO_KEEP & seu$timepoint %in% TP_USE & !is.na(seu$celltype_use))
seu_use <- subset(seu, cells = rownames(seu@meta.data)[cells_keep])

# group filter: require MIN_CELLS_PER_GROUP per (celltype,timepoint,genotype)
md <- as.data.table(seu_use@meta.data, keep.rownames = "cell")
gcnt <- md[, .N, by=.(celltype_use, timepoint, genotype)]
good_groups <- gcnt[N >= MIN_CELLS_PER_GROUP]
md <- md[good_groups, on=.(celltype_use, timepoint, genotype), nomatch=0L]
seu_use <- subset(seu_use, cells = md$cell)

md <- as.data.table(seu_use@meta.data, keep.rownames = "cell")

# -------------------------- expression matrix --------------------------------
Seurat::DefaultAssay(seu_use) <- ASSAY_USE
expr <- Seurat::GetAssayData(seu_use, slot = SLOT_USE)
if (!inherits(expr, "dgCMatrix")) expr <- as(expr, "dgCMatrix")


rc_dna <- function(dna) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(dna)))
}

# ---------------------------- target models ----------------------------------
get_patterns_miRNA_canonical <- function(seq_dna) {
  if (is.na(seq_dna) || nchar(seq_dna) < 8) return(character(0))
  s7 <- substr(seq_dna, 2, 8)  # 2–8
  s6 <- substr(seq_dna, 2, 7)  # 2–7
  c(
    "7mer-m8" = rc_dna(s7),
    "7mer-1A" = paste0("A", rc_dna(s6)),
    "8mer-1A" = paste0("A", rc_dna(s7))
  )
}

get_patterns_piRNA_extended <- function(seq_dna) {
  if (is.na(seq_dna) || nchar(seq_dna) < 12) return(character(0))
  s7  <- substr(seq_dna, 2, 8)
  s10 <- substr(seq_dna, 2, 11)
  s11 <- substr(seq_dna, 2, 12)
  c(
    "seed_2_8" = rc_dna(s7),
    "ext_2_11" = rc_dna(s10),
    "ext_2_12" = rc_dna(s11)
  )
}

get_patterns_offset_kmer <- function(seq_dna, offset = 1, k = 7) {
  if (is.na(seq_dna) || nchar(seq_dna) < (offset + k - 1)) return(character(0))
  s <- substr(seq_dna, offset, offset + k - 1)
  setNames(rc_dna(s), paste0(k, "mer_off", offset))
}

model_def <- list(
  miRNA_canonical = function(seq_dna) {
    pats <- get_patterns_miRNA_canonical(seq_dna)
    w <- c("7mer-m8"=2, "7mer-1A"=1, "8mer-1A"=3)
    list(patterns=pats, weights=w)
  },
  piRNA_extended = function(seq_dna) {
    pats <- get_patterns_piRNA_extended(seq_dna)
    w <- c("seed_2_8"=1, "ext_2_11"=3, "ext_2_12"=4)
    list(patterns=pats, weights=w)
  },
  off1_7mer = function(seq_dna) {
    pats <- get_patterns_offset_kmer(seq_dna, offset=1, k=7)
    w <- setNames(2, names(pats))
    list(patterns=pats, weights=w)
  },
  off2_7mer = function(seq_dna) {
    pats <- get_patterns_offset_kmer(seq_dna, offset=2, k=7)
    w <- setNames(2, names(pats))
    list(patterns=pats, weights=w)
  },
  off3_7mer = function(seq_dna) {
    pats <- get_patterns_offset_kmer(seq_dna, offset=3, k=7)
    w <- setNames(2, names(pats))
    list(patterns=pats, weights=w)
  }
)

# ---------------------- fast scanning with PDict ------------------------------
count_weighted_hits <- function(subject_seqs, patterns_named, weights_named, fixed = TRUE) {
  if (length(patterns_named) == 0) return(setNames(numeric(0), character(0)))
  
  nm <- intersect(names(patterns_named), names(weights_named))
  patterns_named <- patterns_named[nm]
  weights_named  <- weights_named[nm]
  if (length(patterns_named) == 0) return(setNames(numeric(0), character(0)))
  
  subj_names <- names(subject_seqs)
  if (is.null(subj_names)) subj_names <- as.character(seq_along(subject_seqs))
  
  score <- setNames(numeric(length(subject_seqs)), subj_names)
  
  wlen <- nchar(unname(patterns_named))
  idx_by_len <- split(seq_along(patterns_named), wlen)
  
  ns <- length(subject_seqs)
  
  for (idx in idx_by_len) {
    pats <- Biostrings::DNAStringSet(unname(patterns_named[idx]))
    names(pats) <- names(patterns_named)[idx]
    np <- length(pats)
    
    cnt <- tryCatch(
      {
        pd <- Biostrings::PDict(pats)
        Biostrings::vcountPDict(pd, subject_seqs, fixed = fixed)
      },
      error = function(e) {
        # vapply returns subjects × patterns when np > 1
        vapply(seq_along(pats), function(i) {
          Biostrings::vcountPattern(pats[[i]], subject_seqs, fixed = fixed)
        }, FUN.VALUE = integer(ns))
      }
    )
    
    # ---- normalize cnt to patterns × subjects ----
    if (is.null(dim(cnt))) {
      # single pattern: vector of length ns
      cnt <- matrix(as.integer(cnt), nrow = 1, ncol = ns)
    } else {
      cnt <- as.matrix(cnt)
      # if subjects × patterns, transpose
      if (nrow(cnt) == ns && ncol(cnt) == np) cnt <- t(cnt)
    }
    
    # final sanity check
    if (nrow(cnt) != np || ncol(cnt) != ns) {
      stop(sprintf(
        "Unexpected count matrix shape: %d×%d; expected %d×%d (patterns×subjects).",
        nrow(cnt), ncol(cnt), np, ns
      ))
    }
    
    rownames(cnt) <- names(pats)
    colnames(cnt) <- subj_names
    
    w2 <- as.numeric(weights_named[rownames(cnt)])
    w2[is.na(w2)] <- 0
    
    score <- score + as.numeric(matrix(w2, nrow = 1) %*% cnt)
  }
  
  score
}



scan_targets_genelevel <- function(seq_dna,
                                   utr_seqs, cds_seqs, tx2gene_named,
                                   patterns, weights,
                                   use_cds = TRUE,
                                   cds_mult = 0.5,
                                   min_gene_score = 3,
                                   top_n_genes = 400) {
  if (is.na(seq_dna) || length(patterns) == 0) return(data.table(gene=character(0), weight=numeric(0)))
  
  utr_score <- count_weighted_hits(utr_seqs, patterns, weights, fixed=FALSE)
  tx_score <- utr_score
  
  if (isTRUE(use_cds)) {
    cds_score <- count_weighted_hits(cds_seqs, patterns, weights, fixed=FALSE)
    # align / add
    tx_all <- union(names(tx_score), names(cds_score))
    out <- setNames(numeric(length(tx_all)), tx_all)
    out[names(tx_score)] <- out[names(tx_score)] + tx_score
    out[names(cds_score)] <- out[names(cds_score)] + cds_mult * cds_score
    tx_score <- out
  }
  
  tx_keep <- names(tx_score)[tx_score >= min_gene_score]
  if (length(tx_keep) == 0) return(data.table(gene=character(0), weight=numeric(0)))
  
  g <- as.character(tx2gene_named[tx_keep])
  ok <- !is.na(g) & nzchar(g)
  if (!any(ok)) return(data.table(gene=character(0), weight=numeric(0)))
  
  # sum transcript scores per gene
  gene_score <- tapply(tx_score[tx_keep[ok]], g[ok], sum)
  gene_score <- sort(gene_score, decreasing = TRUE)
  gene_score <- gene_score[seq_len(min(length(gene_score), top_n_genes))]
  
  data.table(gene = names(gene_score), weight = as.numeric(gene_score))
}

# ----------------------- build / load target cache ----------------------------
srnas_to_run <- unique(bulk_dir_tbl$sRNA)
srna_seq_tbl <- data.table(sRNA = srnas_to_run)
srna_seq_tbl[, seq_dna := vapply(sRNA, extract_seq_from_id, character(1))]
#if (file.exists(TARGET_CACHE_RDS)) file.remove(TARGET_CACHE_RDS)

if (file.exists(TARGET_CACHE_RDS)) {
  message("Loading target cache: ", TARGET_CACHE_RDS)
  target_cache <- readRDS(TARGET_CACHE_RDS)
} else {
  message("Building targets (can be slow). Will save cache to: ", TARGET_CACHE_RDS)
  
  target_cache <- list()  # names: "<sRNA>__<model>" ; value: data.table(gene, weight)
  for (sid in srna_seq_tbl$sRNA) {
    seq_dna <- srna_seq_tbl[sRNA == sid, seq_dna][1]
    if (is.na(seq_dna)) next
    
    for (m in TARGET_MODELS) {
      def <- model_def[[m]](seq_dna)
      if (length(def$patterns) == 0) next
      
      tg <- scan_targets_genelevel(
        seq_dna = seq_dna,
        utr_seqs = utr_tx_seq,
        cds_seqs = cds_tx_seq,
        tx2gene_named = tx2gene_map,
        patterns = def$patterns,
        weights = def$weights,
        use_cds = USE_CDS_IN_SCAN,
        cds_mult = CDS_WEIGHT_MULT,
        min_gene_score = MIN_GENE_SCORE,
        top_n_genes = TOP_N_GENES
      )
      
      # keep only genes present in scRNA matrix
      tg <- tg[gene %in% rownames(expr)]
      if (nrow(tg) == 0) next
      
      key <- paste0(sid, "__", m)
      target_cache[[key]] <- tg
    }
  }
  saveRDS(target_cache, TARGET_CACHE_RDS)
}

if (length(target_cache) == 0) stop("No targets in cache. Check sequence parsing and scan thresholds.")

#saveRDS(target_cache, "G:/PhD_final/tables/srna_activity_sc/targets_cache_copy.rds")
#TARGET_CACHE_RDS <- file.path(OUT_DIR, "targets_cache.rds")
# ------------------ expression bins for control matching ----------------------
# Use global average expression to bin genes
# Use global average expression to bin genes
gene_means <- Matrix::rowMeans(expr)

# Ensure gene IDs are present as names (robust even if Matrix drops them)
if (is.null(names(gene_means))) names(gene_means) <- rownames(expr)

gene_means <- gene_means[is.finite(gene_means)]
gene_means <- gene_means[intersect(names(gene_means), rownames(expr))]

qs <- quantile(gene_means, probs = seq(0, 1, length.out = N_BINS_EXPR + 1), na.rm = TRUE)
qs <- unique(qs)
if (length(qs) < 5) stop("Expression binning failed (too few unique quantiles).")

gene_bin <- cut(gene_means, breaks = qs, include.lowest = TRUE, labels = FALSE)

# Fix: restore names so split() has gene IDs to split
names(gene_bin) <- names(gene_means)

gene_bin <- gene_bin[!is.na(gene_bin)]
genes_by_bin <- split(names(gene_bin), gene_bin)


seed_from_string <- function(s) {
  x <- utf8ToInt(s)
  as.integer((sum(x) + 131 * length(x)) %% .Machine$integer.max)
}

pick_controls_matched <- function(target_genes, ctrl_per_target = 20L, seed = 1L) {
  tg <- unique(target_genes)
  tg <- tg[tg %in% names(gene_bin)]
  if (length(tg) == 0) return(character(0))
  
  set.seed(seed)
  ctrl <- character(0)
  for (g in tg) {
    b <- gene_bin[[g]]
    pool <- setdiff(genes_by_bin[[as.character(b)]], tg)
    if (length(pool) == 0) next
    take <- min(length(pool), ctrl_per_target)
    ctrl <- c(ctrl, sample(pool, size = take, replace = FALSE))
  }
  unique(ctrl)
}

# ---------------------- activity scoring (per cell) ---------------------------
# Score = (weighted mean target expr) - (mean matched-control expr)
# Activity (repression) = -Score, so that “targets DOWN” => activity increases.
weighted_mean_expr <- function(expr_mat, genes, weights = NULL) {
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr_mat)))
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr_mat)))
  
  sub <- expr_mat[genes, , drop = FALSE]
  if (is.null(weights)) {
    return(Matrix::colMeans(sub))
  }
  w <- weights[match(genes, names(weights))]
  w[is.na(w)] <- 0
  if (sum(w) <= 0) return(Matrix::colMeans(sub))
  w <- w / sum(w)
  
  # t(sub) %*% w  -> vector per cell
  as.numeric(Matrix::t(sub) %*% w)
}

score_srna_activity <- function(expr_mat, tg_dt, ctrl_per_target=20L, min_targets=15L, seed=1L) {
  tg <- unique(tg_dt$gene)
  if (length(tg) < min_targets) return(rep(NA_real_, ncol(expr_mat)))
  
  w <- tg_dt$weight
  names(w) <- tg_dt$gene
  w <- w[w > 0]
  
  ctrl <- pick_controls_matched(tg, ctrl_per_target = ctrl_per_target, seed = seed)
  if (length(ctrl) < min_targets) return(rep(NA_real_, ncol(expr_mat)))
  
  s_tg   <- weighted_mean_expr(expr_mat, tg, weights = w)
  s_ctrl <- weighted_mean_expr(expr_mat, ctrl, weights = NULL)
  
  score <- s_tg - s_ctrl
  activity <- -score
  activity
}

# ---------------------- permutation test helper -------------------------------
perm_test_delta_mean <- function(x, g01, nperm=2000L, seed=1L, alternative=c("two.sided","greater","less")) {
  alternative <- match.arg(alternative)
  ok <- is.finite(x) & !is.na(g01)
  x <- x[ok]
  g01 <- g01[ok]
  if (length(unique(g01)) != 2) return(list(p=NA_real_, obs=NA_real_))
  
  obs <- mean(x[g01 == 1]) - mean(x[g01 == 0])
  
  set.seed(seed)
  n <- length(x)
  g <- as.integer(g01)
  null <- numeric(nperm)
  for (i in seq_len(nperm)) {
    gp <- sample(g, size = n, replace = FALSE)
    null[i] <- mean(x[gp == 1]) - mean(x[gp == 0])
  }
  
  if (alternative == "two.sided") {
    p <- (1 + sum(abs(null) >= abs(obs))) / (nperm + 1)
  } else if (alternative == "greater") {
    p <- (1 + sum(null >= obs)) / (nperm + 1)
  } else {
    p <- (1 + sum(null <= obs)) / (nperm + 1)
  }
  
  list(p=p, obs=obs)
}

# ------------------- expected direction from bulk table -----------------------
if (!("expected_target_change" %in% names(bulk_dir_tbl))) {
  if (!("bulk_log2FC_3dpa" %in% names(bulk_dir_tbl))) stop("bulk_dir_tbl needs expected_target_change or bulk_log2FC_3dpa")
  bulk_dir_tbl[, expected_target_change := ifelse(bulk_log2FC_3dpa > 0, "DOWN_in_ELAC", "UP_in_ELAC")]
}
# expected activity delta sign: DOWN_in_ELAC -> activity higher in ELAC (+)
bulk_dir_tbl[, expected_activity_sign := ifelse(expected_target_change == "DOWN_in_ELAC", +1, -1)]

# ------------------- run: per sRNA × model activity summary -------------------
# Prepare strata ids
md[, strata := paste(celltype_use, timepoint, sep="||")]
md[, geno01 := ifelse(genotype == "ELAC", 1L, 0L)]

# Keep only strata with both genotypes present and enough cells
strata_ok <- md[, .(n_WT = sum(geno01 == 0), n_ELAC = sum(geno01 == 1)), by=strata]
strata_ok <- strata_ok[n_WT >= MIN_CELLS_PER_GROUP & n_ELAC >= MIN_CELLS_PER_GROUP]
md <- md[strata %in% strata_ok$strata]

if (nrow(strata_ok) == 0) stop("No strata pass MIN_CELLS_PER_GROUP for both WT and ELAC.")

# Main results collector
res_list <- list()

keys <- names(target_cache)
message("Scoring ", length(keys), " (sRNA,model) target sets.")


for (key in keys) {
  sid   <- sub("__.*$", "", key)
  model <- sub("^.*__", "", key)
  
  tg_dt <- target_cache[[key]]
  if (is.null(tg_dt) || nrow(tg_dt) < MIN_TARGETS_FOR_SCORE) next
  
  seed_key <- seed_from_string(key)
  
  # score activity
  activity <- score_srna_activity(
    expr_mat = expr,
    tg_dt = tg_dt,
    ctrl_per_target = CTRL_PER_TARGET,
    min_targets = MIN_TARGETS_FOR_SCORE,
    seed = seed_key
  )
  if (is.null(activity) || all(is.na(activity))) next
  
  # attach to metadata (keep only cells that exist in expr)
  dt <- md[, .(cell, celltype_use, timepoint, genotype, geno01, strata)]
  dt <- dt[cell %in% colnames(expr)]
  
  if (!is.null(names(activity))) {
    dt[, activity := activity[cell]]
  } else {
    dt[, activity := activity[match(cell, colnames(expr))]]
  }
  
  # summarize means per group
  sm <- dt[, .(
    n = .N,
    mean_activity = mean(activity, na.rm = TRUE),
    sd_activity = sd(activity, na.rm = TRUE)
  ), by = .(celltype_use, timepoint, genotype, strata)]
  
  # wide tables
  smw  <- data.table::dcast(sm,  celltype_use + timepoint + strata ~ genotype, value.var = "mean_activity")
  smn  <- data.table::dcast(sm,  celltype_use + timepoint + strata ~ genotype, value.var = "n")
  smsd <- data.table::dcast(sm,  celltype_use + timepoint + strata ~ genotype, value.var = "sd_activity")
  
  # need WT and ELAC to compute delta
  if (!all(c("WT", "ELAC") %in% names(smw))) next
  if (!all(c("WT", "ELAC") %in% names(smn))) next
  if (!all(c("WT", "ELAC") %in% names(smsd))) next
  
  smw[, delta_ELAC_minus_WT := ELAC - WT]
  smw[, n_WT   := smn[["WT"]]]
  smw[, n_ELAC := smn[["ELAC"]]]
  
  # effect size (Cohen's d) using pooled SD
  pooled_sd <- sqrt(
    ((smw$n_WT - 1) * (smsd[["WT"]]^2) + (smw$n_ELAC - 1) * (smsd[["ELAC"]]^2)) /
      pmax(smw$n_WT + smw$n_ELAC - 2, 1)
  )
  pooled_sd <- pmax(pooled_sd, 1e-8)
  smw[, cohen_d := delta_ELAC_minus_WT / pooled_sd]
  
  # expected direction info from bulk (3 dpa)
  exp <- bulk_dir_tbl[sRNA == sid]
  if (nrow(exp) == 0) next
  
  exp_sign <- exp$expected_activity_sign[1]
  smw[, expected_activity_sign := exp_sign]
  smw[, expected_activity_change := ifelse(exp_sign > 0, "DOWN_in_ELAC", "UP_in_ELAC")]
  smw[, bulk_log2FC_3dpa := exp$bulk_log2FC_3dpa[1]]
  smw[, expected_target_change := exp$expected_target_change[1]]
  
  # expected-direction filter at stratum level
  smw[, delta_expected := delta_ELAC_minus_WT * expected_activity_sign]
  smw[, ok_dir := is.finite(delta_expected) & (delta_expected > 0)]
  
  # permutation p-values
  smw[, perm_p_two_sided := NA_real_]
  smw[, perm_p_expected  := NA_real_]
  
  if (isTRUE(DO_PERMUTATION)) {
    
    strata_to_test <- unique(smw[ok_dir == TRUE, strata])
    
    for (st in strata_to_test) {
      sub <- dt[strata == st & is.finite(activity)]
      
      # require both groups to have enough cells
      n0 <- sub[, sum(geno01 == 0L)]
      n1 <- sub[, sum(geno01 == 1L)]
      if (n0 < MIN_CELLS_PER_GROUP || n1 < MIN_CELLS_PER_GROUP) next
      
      # two-sided
      pt2 <- perm_test_delta_mean(
        x = sub$activity,
        g01 = sub$geno01,
        nperm = N_PERM,
        seed = seed_key + seed_from_string(st),
        alternative = "two.sided"
      )
      
      # directional (expected sign)
      alt <- if (exp_sign > 0) "greater" else "less"
      ptd <- perm_test_delta_mean(
        x = sub$activity,
        g01 = sub$geno01,
        nperm = N_PERM,
        seed = seed_key + 7L + seed_from_string(st),
        alternative = alt
      )
      
      smw[strata == st, perm_p_two_sided := pt2$p]
      smw[strata == st, perm_p_expected  := ptd$p]
    }
  }
  
  # annotate and store
  smw[, sRNA := sid]
  smw[, model := model]
  smw[, n_targets := nrow(tg_dt)]
  
  res_list[[length(res_list) + 1L]] <- smw
}

res <- rbindlist(res_list, fill = TRUE)
if (nrow(res) == 0) stop("No results produced. Check target cache + scoring thresholds.")

# multiple testing correction within each (sRNA,model) over strata
res[, perm_p_expected_BH := p.adjust(perm_p_expected, method = "BH"), by=.(sRNA, model)]
res[, perm_p_two_sided_BH := p.adjust(perm_p_two_sided, method = "BH"), by=.(sRNA, model)]

# prioritization score that respects expected direction (only meaningful where perm_p_expected computed)
res[, score_expected := expected_activity_sign * delta_ELAC_minus_WT * (-log10(perm_p_expected_BH + 1e-300))]

# write outputs
#fwrite(res, file.path(OUT_DIR, "srna_activity_moduleScore_controlMatched_72_only_correct_delta.tsv"), sep = "\t")
#save(res,file="G:/PhD_final/tables/srna_activity_moduleScore_controlMatched_72_only_correct_delta.RData")
#save(activity,file="G:/PhD_final/tables/srna_activity_72_only_correct_delta.RData")
# “best hit” per (sRNA,model) using BH directional p, then |delta|
best <- res[order(perm_p_expected_BH, -abs(delta_ELAC_minus_WT))]
best <- best[, .SD[1], by=.(sRNA, model)]
#fwrite(best, file.path(OUT_DIR, "srna_activity_bestHit.tsv"), sep = "\t")

best$consistensy <- ifelse(((best$ELAC<best$GFP) & (best$ELAC<best$WT)) |
                             ((best$ELAC>best$GFP) & (best$ELAC>best$WT)),
                           TRUE,FALSE)
best[best$consistensy==TRUE,]
#fwrite(best, file.path(OUT_DIR, "srna_activity_bestHit.tsv"), sep = "\t")
best[best$sRNA=="GCATCGGTGGTTCAGTGGTAGAATGCTCGCCT 5'-tiRNA-Gly-GCC",]
load("G:/PhD_final/tables/srna_activity_72_only_correct_delta.RData")#activity

library(scCustomize)
result_obj


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


################################################################################
# Rescaled plot_sncRNA_activity_umap
?scale_fill_gradient2
plot_sncRNA_activity_umap_rescaled <- function(
    seu,
    activity,
    sid,
    model,
    focus_pops = c("Phagocyte (broad)", "PGRN⁺ parenchymal cell"),
    detailed_cols,
    reduction = "umap.d33.nn100.md0.3",
    pop_col = "final_population",
    cond_col = "condition_correct",
    facet_map = NULL,
    facet_levels = NULL,
    title = NULL,
    rel_widths = c(1, 2),
    
    # palette (will be pastelized by default)
    low_col  = "#A6CEE3",
    mid_col  = "#FFFFFF",
    high_col = "#FB9A99",
    pastelize = TRUE,
    pastel_strength = 0.2,   # 0..1, higher = lighter
    
    # point sizes
    pt_bg = 0.10,
    pt_pop = 0.25,
    pt_act = 0.55,
    
    # alphas and background
    alpha_bg = 0.45,
    alpha_act = 0.95,
    bg_col = "grey88",
    
    # legends + text
    legend_dot_size = 3,
    strip_text_size = 11,
    legend_text_size  = 12,
    legend_title_size = 13,
    title_size        = 14,
    title_gap_lines   = 6,
    
    # activity visualization controls
    activity_mode = c("raw", "delta_vs_ref"),
    ref_condition = "WT72",
    delta_within = c("pop", "global"),
    act_clip_q = c(0.05, 0.95),
    act_limits = NULL,
    diverging_symmetric = TRUE,
    midpoint = 0,
    
    # transform (new ggplot2 arg name is transform, not trans)
    act_transform = c("auto", "identity", "modulus"),
    modulus_p = 0.7,
    
    # drawing style
    use_fill = FALSE,         # default FALSE to avoid dark outlines from shape 21
    outline_col = NA,         # NA removes border for shape 21
    outline_stroke = 0,
    order_by_abs = TRUE
) {
  # ---- helpers ----
  .mix_with_white <- function(col, w = 0.55) {
    if (is.na(col) || !nzchar(col)) return(col)
    rgb <- grDevices::col2rgb(col) / 255
    rgb2 <- rgb * (1 - w) + 1 * w
    grDevices::rgb(rgb2[1], rgb2[2], rgb2[3])
  }
  
  # ---- checks ----
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(!is.null(names(activity)))
  stopifnot(is.character(sid), length(sid) == 1)
  stopifnot(is.character(model), length(model) == 1)
  stopifnot(is.character(reduction), length(reduction) == 1)
  stopifnot(is.character(pop_col), length(pop_col) == 1)
  stopifnot(is.character(cond_col), length(cond_col) == 1)
  
  if (missing(detailed_cols) || is.null(detailed_cols)) {
    stop("Provide 'detailed_cols' (named vector: population -> color).")
  }
  if (!all(focus_pops %in% names(detailed_cols))) {
    stop("Some focus_pops are missing from names(detailed_cols): ",
         paste(setdiff(focus_pops, names(detailed_cols)), collapse = ", "))
  }
  
  activity_mode <- match.arg(activity_mode)
  delta_within  <- match.arg(delta_within)
  act_transform <- match.arg(act_transform)
  
  # ---- pastelize palette (reduces “dark” look on dense UMAPs) ----
  if (isTRUE(pastelize)) {
    w <- max(0, min(1, pastel_strength))
    low_col  <- .mix_with_white(low_col,  w = w)
    mid_col  <- .mix_with_white(mid_col,  w = w / 2)   # keep midpoint close to white
    high_col <- .mix_with_white(high_col, w = w)
  }
  
  # ---- activity column name ----
  act_col <- paste0("act__", make.names(paste(sid, model, sep = "__")))
  
  # ---- align activity to Seurat cells and store ----
  cells <- Seurat::Cells(seu)
  act_vec <- activity[cells]
  names(act_vec) <- cells
  
  n_match <- sum(is.finite(act_vec))
  if (n_match == 0) {
    stop("No overlap between names(activity) and Cells(seu). ",
         "Check: head(names(activity)) vs head(Cells(seu)).")
  }
  seu[[act_col]] <- act_vec
  
  # ---- embeddings ----
  if (!reduction %in% names(seu@reductions)) {
    stop("Reduction '", reduction, "' not found. Available: ",
         paste(names(seu@reductions), collapse = ", "))
  }
  um <- Seurat::Embeddings(seu, reduction = reduction)[, 1:2, drop = FALSE]
  colnames(um) <- c("UMAP_1", "UMAP_2")
  
  df <- data.table::as.data.table(um, keep.rownames = "cell")
  
  meta <- Seurat::FetchData(seu, vars = c(pop_col, cond_col, act_col))
  meta$cell <- rownames(meta)
  
  df <- data.table::merge.data.table(
    df,
    data.table::as.data.table(meta),
    by = "cell"
  )
  
  data.table::setnames(df, pop_col, "pop")
  data.table::setnames(df, cond_col, "cond")
  data.table::setnames(df, act_col, "activity_raw")
  
  # ---- facet labels (plotmath strings) ----
  if (is.null(facet_map)) {
    facet_map <- c(
      "ELAC72" = "paste(italic('Smed ELAC2'), ' KD 72 hpa')",
      "GFP72"  = "'GFP mock 72 hpa'",
      "WT72"   = "'WT 72 hpa'"
    )
  }
  
  df[, condition_facet := {
    cc <- as.character(cond)
    out <- unname(facet_map[cc])
    out[is.na(out)] <- shQuote(cc[is.na(out)])
    out
  }]
  
  if (is.null(facet_levels)) {
    base_levels <- unname(facet_map)
    extras <- setdiff(unique(df$condition_facet), base_levels)
    facet_levels <- c(base_levels, extras)
  }
  df[, condition_facet := factor(condition_facet, levels = facet_levels)]
  
  # ---- focus flags ----
  df[, is_focus := pop %in% focus_pops]
  
  # ---- compute activity to plot (raw or delta) ----
  df[, activity_plot := as.numeric(NA)]
  df[is_focus == TRUE, activity_plot := activity_raw]
  
  df_focus <- df[is_focus == TRUE & is.finite(activity_plot)]
  if (nrow(df_focus) == 0) stop("No finite activity values found within focus_pops.")
  
  if (activity_mode == "delta_vs_ref") {
    if (!ref_condition %in% unique(df$cond)) {
      stop("ref_condition '", ref_condition, "' not found in cond. Present: ",
           paste(sort(unique(df$cond)), collapse = ", "))
    }
    
    ref_dt <- df[is_focus == TRUE & cond == ref_condition & is.finite(activity_plot)]
    if (nrow(ref_dt) == 0) {
      stop("No finite focus activity values in ref_condition = '", ref_condition, "'.")
    }
    
    if (delta_within == "pop") {
      base <- ref_dt[, .(baseline = stats::median(activity_plot, na.rm = TRUE)), by = pop]
      df <- data.table::merge.data.table(df, base, by = "pop", all.x = TRUE)
      df[is_focus == TRUE & is.finite(activity_plot), activity_plot := activity_plot - baseline]
      df[, baseline := NULL]
    } else {
      baseline <- stats::median(ref_dt$activity_plot, na.rm = TRUE)
      df[is_focus == TRUE & is.finite(activity_plot), activity_plot := activity_plot - baseline]
    }
    
    df_focus <- df[is_focus == TRUE & is.finite(activity_plot)]
  }
  
  # ---- order so extremes draw on top ----
  if (order_by_abs) df_focus <- df_focus[order(abs(activity_plot))]
  
  # ---- robust limits (clip + squish) ----
  vals <- df_focus$activity_plot
  
  if (!is.null(act_limits)) {
    if (!is.numeric(act_limits) || length(act_limits) != 2) {
      stop("act_limits must be numeric length 2, e.g. c(-0.2, 0.2).")
    }
    lims <- sort(as.numeric(act_limits))
  } else {
    qq <- stats::quantile(vals, probs = act_clip_q, na.rm = TRUE, names = FALSE)
    lims <- sort(as.numeric(qq))
  }
  
  if (diverging_symmetric) {
    max_abs <- max(abs(lims - midpoint))
    lims <- midpoint + c(-max_abs, max_abs)
  }
  
  # ---- pick transform object ----
  act_transform_obj <- switch(
    act_transform,
    "identity" = "identity",
    "modulus"  = scales::transform_modulus(p = modulus_p),
    "auto"     = if (activity_mode == "delta_vs_ref") scales::transform_modulus(p = modulus_p) else "identity"
  )
  
  # ---- right panel background data ----
  df_nonfocus <- df[is_focus == FALSE]
  df_focus_na <- df[is_focus == TRUE & !is.finite(activity_plot)]
  
  # ---- left: highlight focus pops ----
  p_left <- ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(color = bg_col, size = pt_bg, alpha = alpha_bg) +
    ggplot2::geom_point(
      data = df[df$pop %in% focus_pops],
      ggplot2::aes(color = pop),
      size = pt_pop,
      alpha = 1
    ) +
    ggplot2::scale_color_manual(
      values = detailed_cols[focus_pops],
      breaks = focus_pops,
      name = "Population"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = legend_dot_size, alpha = 1))
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text  = ggplot2::element_text(size = legend_text_size),
      text = ggplot2::element_text(size = 11)
    )
  
  # ---- right: non-focus grey, focus colored by activity_plot ----
  scale_name <- if (activity_mode == "delta_vs_ref") "Δ Activity" else "Activity"
  
  if (use_fill) {
    p_right <- ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2)) +
      ggplot2::geom_point(data = df_nonfocus, color = bg_col, size = pt_bg, alpha = alpha_bg) +
      ggplot2::geom_point(data = df_focus_na, color = bg_col, size = pt_bg, alpha = alpha_bg) +
      ggplot2::geom_point(
        data = df_focus,
        ggplot2::aes(fill = activity_plot),
        shape = 21,
        color = outline_col,
        stroke = outline_stroke,
        size = pt_act,
        alpha = alpha_act
      ) +
      ggplot2::scale_fill_gradient2(
        name = scale_name,
        low = low_col, mid = mid_col, high = high_col,
        midpoint = midpoint,
        limits = lims,
        oob = scales::squish,
        transform = act_transform_obj
      ) +
      ggplot2::facet_wrap(~ condition_facet, nrow = 1, labeller = ggplot2::label_parsed) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 18, b = 0, l = 0),  # more space on the right
        legend.box.margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 6),
        legend.title = ggplot2::element_text(size = legend_title_size),
        legend.text  = ggplot2::element_text(size = legend_text_size),
        strip.text   = ggplot2::element_text(size = strip_text_size)
      )
      # ggplot2::theme(
      #   legend.position = "right",
      #   legend.title = ggplot2::element_text(size = legend_title_size),
      #   legend.text  = ggplot2::element_text(size = legend_text_size),
      #   strip.text   = ggplot2::element_text(size = strip_text_size)
      # )
  } else {
    p_right <- ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2)) +
      ggplot2::geom_point(data = df_nonfocus, color = bg_col, size = pt_bg, alpha = alpha_bg) +
      ggplot2::geom_point(data = df_focus_na, color = bg_col, size = pt_bg, alpha = alpha_bg) +
      ggplot2::geom_point(
        data = df_focus,
        ggplot2::aes(color = activity_plot),
        size = pt_act,
        alpha = alpha_act
      ) +
      ggplot2::scale_color_gradient2(
        name = scale_name,
        low = low_col, mid = mid_col, high = high_col,
        midpoint = midpoint,
        limits = lims,
        oob = scales::squish,
        transform = act_transform_obj
      ) +
      ggplot2::facet_wrap(~ condition_facet, nrow = 1, labeller = ggplot2::label_parsed) +
      ggplot2::theme_void() +
      # ggplot2::theme(
      #   legend.position = "right",
      #   legend.title = ggplot2::element_text(size = legend_title_size),
      #   legend.text  = ggplot2::element_text(size = legend_text_size),
      #   strip.text   = ggplot2::element_text(size = strip_text_size)
      # )
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 18, b = 0, l = 0),  # more space on the right
        legend.box.margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 6),
        legend.title = ggplot2::element_text(size = legend_title_size),
        legend.text  = ggplot2::element_text(size = legend_text_size),
        strip.text   = ggplot2::element_text(size = strip_text_size)
      )
  }
  
  # ---- combine + title ----
  core <- cowplot::plot_grid(
    p_left, p_right,
    ncol = 2,
    rel_widths = rel_widths,
    align = "h"
  )
  
  if (is.null(title)) title <- paste0(sid, " | ", model)
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(
      title,
      x = 0.05, y = 1,
      hjust = 0, vjust = 1.25,
      fontface = "bold",
      size = title_size
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 0, r = 0, b = title_gap_lines, l = 0)
    )
  
  final <- cowplot::plot_grid(
    title_grob,
    core,
    ncol = 1,
    rel_heights = c(0.10, 1)
  )
  
  list(
    plot = final,
    p_left = p_left,
    p_right = p_right,
    act_col = act_col,
    df = df,
    df_focus = df_focus,
    seu = seu,
    activity_limits = lims,
    activity_mode = activity_mode,
    ref_condition = ref_condition,
    act_transform = act_transform_obj
  )
}

out1 <- plot_sncRNA_activity_umap_rescaled(
  seu = seu_use,
  activity = activity,
  sid = "GCATCGGTGGTTCAGTGGTAGAATGCTCGCCT 5'-tiRNA-Gly-GCC",
  model = "miRNA_canonical",
  focus_pops = c("PGRN⁺ parenchymal cell", "Phagocyte (broad)"),
  detailed_cols = DETAILED_COLS,
  title = "5'-tiRNA-Gly-GCC (GCATCGGTGGTTCAGTGGTAGAATGCTCGCCT)",
  activity_mode = "delta_vs_ref",
  ref_condition = "WT72",
  delta_within = "pop",
  act_clip_q = c(0.10, 0.90),
  low_col = "blue",
  high_col = "red",
  mid_col = "white"
)
#out1_new$plot



out2 <- plot_sncRNA_activity_umap_rescaled(
  seu = seu_use,
  activity = activity,
  sid = "TCTTTGGTTTTCTAGC sme-miR-9a-5p",
  model = "off1_7mer",
  focus_pops = c("Glutamatergic neuron"),
  detailed_cols = DETAILED_COLS,
  title = "miR-9a-5p (TCTTTGGTTTTCTAGC)",
  activity_mode = "delta_vs_ref",
  ref_condition = "WT72",
  delta_within = "pop",
  act_clip_q = c(0.10, 0.90),
  low_col = "blue",
  high_col = "red",
  mid_col = "white"
)


out4 <- plot_sncRNA_activity_umap_rescaled(
  seu = seu_use,
  activity = activity,
  sid = "CACATGACATGTATACTCTACAAACGCAC piRNA",
  model = "piRNA_extended",
  focus_pops = c("ζ-neoblast (epidermal-fated)"),
  detailed_cols = DETAILED_COLS,
  title = "piRNA (CACATGACATGTATACTCTACAAACGCAC)",
  activity_mode = "delta_vs_ref",
  ref_condition = "WT72",
  delta_within = "pop",
  act_clip_q = c(0.10, 0.90),
  low_col = "blue",
  high_col = "red",
  mid_col = "white"
)

out3 <- plot_sncRNA_activity_umap_rescaled(
  seu = seu_use,
  activity = activity,
  sid = "ACCACTGACCGAGCATATCC sme-miR-190a-3p",
  model = "off1_7mer",
  focus_pops = c("Late epidermal progenitor"),
  detailed_cols = DETAILED_COLS,
  title = "miR-190a-3p (ACCACTGACCGAGCATATCC)",
  activity_mode = "delta_vs_ref",
  ref_condition = "WT72",
  delta_within = "pop",
  act_clip_q = c(0.10, 0.90),
  low_col = "blue",
  high_col = "red",
  mid_col = "white"
)


library(ggpubr)
ggarrange(out1$plot,
          out2$plot,
          out3$plot,
          out4$plot,
          ncol = 1,
          labels = c("A","B","C","D"),
          common.legend = TRUE)

