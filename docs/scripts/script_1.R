
#   Clean, reproducible pipeline to:
#   1) read Bowtie2-mapped small-RNA tables from multiple references (RNAcentral,
#      planarian tRNA set, rRNA set, transcriptome, genome),
#   2) assign an RNA class (miRNA, tRF, rRF, etc.) using consistent rules,
#   3) compute a *consensus* annotation per unique sequence (resolve multi-hits),
#   4) build a sample-by-feature count matrix and CPM-filter it for DGE.


# ----------------------------
# 0) Packages 
# ----------------------------
install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!bioc) {
        install.packages(p, repos = "https://cloud.r-project.org")
      } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(p, update = FALSE, ask = FALSE)
      }
    }
  }
}

cran_pkgs <- c("data.table", "dplyr", "stringr", "purrr", "tidyr", "forcats", "ggplot2", "ggrepel", "scales")
bioc_pkgs <- c("GenomicAlignments", "Biostrings", "edgeR")

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

# ----------------------------
# 1) Configuration
# ----------------------------
cfg <- list(
  # Input directories
  dir_rnacentral  = "F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_RNAcentral/mapped_seq_with_strand_new/",
  dir_trna_map    = "E:/Illumina/PARN_ELAC2_silencing/smallRNA/smallRNAwithAdapters/miRNA/bowtie2_mapped_SM_tRNA_seq_with_strand/",
  dir_rrna_map    = "F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_rRNA_and_genome/rRNA_seq_with_strand/",
  dir_transcript  = "F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_transcriptome/mapped_seq_with_strand/",
  dir_genome      = "F:/PARN_ELAC_silencing/smallRNA/calculated_data_bowtie2_end_to_end/mapped_to_genome/mapped_seq_with_strand/",
  
  # Output
  out_dir         = "F:/PARN_ELAC_silencing/smallRNA/clean_pipeline_out/",
  
  # Filtering
  keep_strand     = 0L,   # keep only strand == 0 
  max_total_mm    = 3L,   # keep reads with (XM + non-match cigar ops) <= 3
  
  # CPM filtering for DGE
  cpm_threshold   = 10,
  cpm_min_samples = 2L
)

if (!dir.exists(cfg$out_dir)) dir.create(cfg$out_dir, recursive = TRUE)

# ----------------------------
# 2) Helpers
# ----------------------------


revcomp_chr <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

# Read a mapping table 
read_map_tbl <- function(path) {
  dt <- data.table::fread(
    path,
    header = FALSE,
    fill = TRUE,
    showProgress = FALSE
  )
  
  # Keep only the first 9 columns and name them consistently.
  if (ncol(dt) < 9) {
    stop("File has < 9 columns: ", path)
  }
  dt <- dt[, 1:9]
  data.table::setnames(dt, c("read", "strand", "ref", "position", "qual", "cigar", "seq", "XM", "MD"))
  dt
}

parse_xm <- function(xm_chr) {
  out <- suppressWarnings(as.integer(stringr::str_replace(xm_chr, "^XM:i:", "")))
  out[is.na(out)] <- 0L
  out
}

# Compute non-match operations in CIGAR using GenomicAlignments::cigarOpTable.
cigar_nonmatch_ops <- function(cigar_chr) {
  op <- GenomicAlignments::cigarOpTable(cigar_chr)
  
  # Count matches separately (M plus "=" if present); everything else is "non-match ops".
  match_cols <- intersect(colnames(op), c("M", "="))
  nonmatch_cols <- setdiff(colnames(op), match_cols)
  
  nonmatch <- if (length(nonmatch_cols) == 0) rep.int(0L, nrow(op)) else rowSums(op[, nonmatch_cols, drop = FALSE])
  match    <- if (length(match_cols) == 0) rep.int(0L, nrow(op)) else rowSums(op[, match_cols, drop = FALSE])
  
  list(nonmatch = as.integer(nonmatch), match = as.integer(match))
}

# RNA class assignment (single source of truth, used everywhere).
rna_type_from_ref <- function(ref_chr) {
  # Priority is encoded by order of case_when clauses (top wins).
  dplyr::case_when(
    stringr::str_detect(ref_chr, stringr::regex("multiple_hits", ignore_case = TRUE)) ~ "multiple_hits",
    stringr::str_detect(ref_chr, stringr::regex("dd_Smed_v6", ignore_case = TRUE)) ~ "mRNA fragments",
    
    # rRNA-like
    stringr::str_detect(ref_chr, stringr::regex("ribosomal|\\brRNA\\b|ITS1|ITS2|SpacerA|\\b28S\\b|\\b12S\\b|\\b16S\\b|5\\.8S|Schmed_cloneH735c", ignore_case = TRUE)) ~ "rRNA fragments",
    
    # tRNA-like
    stringr::str_detect(ref_chr, stringr::regex("\\btRNA\\b|transfer", ignore_case = TRUE)) ~ "tRNA fragments",
    
    # sno/sn
    stringr::str_detect(ref_chr, stringr::regex("snoRNA|nucleolar", ignore_case = TRUE)) ~ "snoRNA fragments",
    stringr::str_detect(ref_chr, stringr::regex("snRNA|spliceosomal|\\b7SK\\b|nuclear", ignore_case = TRUE)) ~ "snRNA fragments",
    
    # small regulatory
    stringr::str_detect(ref_chr, stringr::regex("piRNA", ignore_case = TRUE)) ~ "piRNA",
    stringr::str_detect(ref_chr, stringr::regex("miRNA|microRNA|\\bmiR\\b|Sme-|sme-lin|sme-let|Sme-Bantam", ignore_case = TRUE)) ~ "miRNA",
    
    # lnc
    stringr::str_detect(ref_chr, stringr::regex("long_non-coding|\\blnc\\b", ignore_case = TRUE)) ~ "lncRNA fragments",
    
    TRUE ~ "other fragments"
  )
}

#  sample are ELAC13S, WT25S, GFP33S, etc.
sample_to_group <- function(sample_id) {
  #  ELAC13S -> ELAC3S
  stringr::str_replace(sample_id, "^(ELAC|WT|GFP|PARN)\\d(3S|5S)$", "\\1\\2")
}

list_basenames <- function(dir_path) {
  list.files(dir_path, full.names = FALSE, all.files = FALSE)
}

common_basenames <- function(dirs) {
  bn <- purrr::map(dirs, list_basenames)
  Reduce(intersect, bn)
}

# ----------------------------
# 3) Per-sample processing
# ----------------------------
process_one_sample <- function(basename_file, cfg) {
  # Read each source
  dt_rc   <- read_map_tbl(file.path(cfg$dir_rnacentral, basename_file))
  dt_trna <- read_map_tbl(file.path(cfg$dir_trna_map, basename_file))
  dt_rrna <- read_map_tbl(file.path(cfg$dir_rrna_map, basename_file))
  dt_tx   <- read_map_tbl(file.path(cfg$dir_transcript, basename_file))
  
  # Annotate tRNA/rRNA sources 
  dt_trna[, RNA_type := "tRNA fragments"]
  dt_rrna[, RNA_type := "rRNA fragments"]
  
  # Annotate RNAcentral
  dt_rc[, RNA_type := rna_type_from_ref(ref)]
  
  # Priority logic 
  # 1) keep tRNA mappings first (remove those reads from RNAcentral set)
  dt_rc2 <- dt_rc[!(read %in% dt_trna$read)]
  dt_all <- data.table::rbindlist(list(dt_rc2, dt_trna), use.names = TRUE, fill = TRUE)
  
  # 2) add rRNA mappings only for reads not already present
  dt_rrna2 <- dt_rrna[!(read %in% dt_all$read)]
  dt_all <- data.table::rbindlist(list(dt_all, dt_rrna2), use.names = TRUE, fill = TRUE)
  
  # 3) add transcriptome mappings (mRNA fragments) only for remaining reads
  dt_tx2 <- dt_tx[!(read %in% dt_all$read)]
  dt_tx2[, RNA_type := "mRNA fragments"]
  dt_all <- data.table::rbindlist(list(dt_all, dt_tx2), use.names = TRUE, fill = TRUE)
  
  # Compute mismatch features
  cig <- cigar_nonmatch_ops(dt_all$cigar)
  dt_all[, cigar_match    := cig$match]
  dt_all[, cigar_nonmatch := cig$nonmatch]
  dt_all[, xm_mm          := parse_xm(XM)]
  dt_all[, all_mm         := xm_mm + cigar_nonmatch]
  
  # Compute oriented sequence 
  dt_all[, orig_seq := ifelse(strand == 0, seq, revcomp_chr(seq))]
  
  # Apply filters 
  dt_all <- dt_all[strand == cfg$keep_strand & all_mm <= cfg$max_total_mm]
  
  dt_all[, .(
    read, strand, ref, position, cigar, seq, orig_seq,
    RNA_type, XM, MD, cigar_match, cigar_nonmatch, xm_mm, all_mm
  )]
}

# ----------------------------
# 4) Consensus annotation per unique sequence
# ----------------------------
#   - build unique (orig_seq, pos_ref) pairs from all samples,
#   - classify each pos_ref,
#   - if a sequence maps to multiple types, prefer:
#       (1) tRNA fragments
#       (2) entries containing "Schmidtea"
#       (3) "other fragments" (only if no tRNA)
#       else "multiple_hits"
build_consensus_annotation <- function(dt_all_samples) {
  dt <- data.table::as.data.table(dt_all_samples)
  
  # pos_ref includes coordinate + ref 
  dt[, pos_ref := paste(position, ref, sep = "_")]
  
  # Unique mapping candidates
  uniq <- unique(dt[, .(orig_seq, pos_ref)])
  
  # Classify based on pos_ref content 
  uniq[, RNA_type := rna_type_from_ref(pos_ref)]
  
  # Determine if any sequence has multiple RNA types
  uniq[, correct_annotation := {
    # All candidates for this orig_seq
    pos_all <- pos_ref
    type_all <- RNA_type
    
    # Preference 1: any tRNA hit
    if (any(type_all == "tRNA fragments")) {
      pos_all[which(type_all == "tRNA fragments")[1]]
    } else if (any(stringr::str_detect(pos_all, stringr::regex("Schmidtea", ignore_case = TRUE)))) {
      pos_all[which(stringr::str_detect(pos_all, stringr::regex("Schmidtea", ignore_case = TRUE)))[1]]
    } else if (any(type_all == "other fragments")) {
      pos_all[which(type_all == "other fragments")[1]]
    } else {
      "multiple_hits"
    }
  }, by = orig_seq]
  
  uniq[, correct_RNA_type := rna_type_from_ref(correct_annotation)]
  
  # One row per orig_seq
  out <- unique(uniq[, .(orig_seq, correct_annotation, correct_RNA_type)])
  out
}

# ----------------------------
# 5) Genome-only add-in 
# ----------------------------
read_genome_only <- function(basename_file, cfg) {
  dt_g <- read_map_tbl(file.path(cfg$dir_genome, basename_file))
  dt_g[, RNA_type := "genome"]
  
  cig <- cigar_nonmatch_ops(dt_g$cigar)
  dt_g[, cigar_match    := cig$match]
  dt_g[, cigar_nonmatch := cig$nonmatch]
  dt_g[, xm_mm          := parse_xm(XM)]
  dt_g[, all_mm         := xm_mm + cigar_nonmatch]
  dt_g[, orig_seq       := ifelse(strand == 0, seq, revcomp_chr(seq))]
  
  dt_g <- dt_g[strand == cfg$keep_strand & all_mm <= cfg$max_total_mm]
  
  dt_g[, .(
    read, strand, ref, position, cigar, seq, orig_seq,
    RNA_type, XM, MD, cigar_match, cigar_nonmatch, xm_mm, all_mm
  )]
}

# ----------------------------
# 6) Count matrix + CPM filtering
# ----------------------------
make_count_matrix <- function(sample_tables) {
  # Build (seq_anno, sample, n) long table and cast wide.
  dt_long <- data.table::rbindlist(
    lapply(names(sample_tables), function(sid) {
      dt <- data.table::as.data.table(sample_tables[[sid]])
      dt[, seq_anno := paste(orig_seq, correct_annotation, sep = " ")]
      dt[, .(n = .N), by = .(seq_anno)][, sample := sid][]
    }),
    use.names = TRUE, fill = TRUE
  )
  
  dt_wide <- data.table::dcast(dt_long, seq_anno ~ sample, value.var = "n", fill = 0)
  
  # Convert to numeric matrix for edgeR::cpm
  mat <- as.matrix(dt_wide[, -1])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- dt_wide$seq_anno
  
  mat
}

cpm_filter <- function(count_mat, cpm_threshold = 10, min_samples = 2L) {
  cpm <- edgeR::cpm(count_mat)
  keep <- rowSums(cpm > cpm_threshold) >= min_samples
  list(cpm = cpm, keep = keep)
}

# ==============================================================================
# 7) Run the pipeline
# ==============================================================================
dirs_needed <- c(cfg$dir_rnacentral, cfg$dir_trna_map, cfg$dir_rrna_map, cfg$dir_transcript)
basenames <- common_basenames(dirs_needed)

if (length(basenames) == 0) {
  stop("No common basenames found across the required directories.")
}

# --- A) Build per-sample tables (no genome)
sample_tbls <- purrr::set_names(vector("list", length(basenames)), basenames)
for (bn in basenames) {
  message("Processing (no genome): ", bn)
  sample_tbls[[bn]] <- process_one_sample(bn, cfg)
}

# Create a stable sample_id from basename 
sample_ids <- stringr::str_replace(basenames, "\\.[^.]*$", "")
names(sample_tbls) <- sample_ids

# Add group labels like ELAC3S / WT5S
for (sid in names(sample_tbls)) {
  dt <- data.table::as.data.table(sample_tbls[[sid]])
  dt[, set := sample_to_group(sid)]
  sample_tbls[[sid]] <- dt
}

# --- B) Build consensus annotation map from all samples
dt_all <- data.table::rbindlist(sample_tbls, use.names = TRUE, fill = TRUE)
anno_map <- build_consensus_annotation(dt_all)

# Save annotation map
saveRDS(anno_map, file.path(cfg$out_dir, "consensus_annotation_map.rds"))
data.table::fwrite(anno_map, file.path(cfg$out_dir, "consensus_annotation_map.csv"))

# --- C) Attach consensus annotation to each sample table
sample_tbls_anno <- lapply(sample_tbls, function(dt) {
  dt <- data.table::as.data.table(dt)
  dt2 <- merge(dt, anno_map, by = "orig_seq", all.x = TRUE)
  
  # If something is NA ...
  dt2[is.na(correct_annotation), correct_annotation := paste(position, ref, sep = "_")]
  dt2[is.na(correct_RNA_type),   correct_RNA_type   := rna_type_from_ref(correct_annotation)]
  
  dt2
})

saveRDS(sample_tbls_anno, file.path(cfg$out_dir, "sample_tables_no_genome_annotated.rds"))

# --- Add genome-only reads (heavy)

basenames_g <- common_basenames(c(cfg$dir_genome, cfg$dir_rnacentral))
basenames_g <- intersect(basenames_g, basenames)  # keep aligned set

genome_tbls <- purrr::set_names(vector("list", length(basenames_g)), basenames_g)
for (bn in basenames_g) {
  message("Processing (genome): ", bn)
  genome_tbls[[bn]] <- read_genome_only(bn, cfg)
}
names(genome_tbls) <- stringr::str_replace(basenames_g, "\\.[^.]*$", "")

# Append genome reads not already assigned in no-genome set (by read ID)
sample_tbls_with_genome <- lapply(names(sample_tbls_anno), function(sid) {
  dt_no <- data.table::as.data.table(sample_tbls_anno[[sid]])
  
  if (!sid %in% names(genome_tbls)) return(dt_no)
  
  dt_g <- data.table::as.data.table(genome_tbls[[sid]])
  dt_g <- dt_g[!(read %in% dt_no$read)]
  dt_g[, set := sample_to_group(sid)]
  dt_g[, correct_annotation := ref]
  dt_g[, correct_RNA_type := "genome"]
  
  data.table::rbindlist(list(dt_no, dt_g), use.names = TRUE, fill = TRUE)
})
names(sample_tbls_with_genome) <- names(sample_tbls_anno)

saveRDS(sample_tbls_with_genome, file.path(cfg$out_dir, "sample_tables_with_genome_annotated.rds"))

# --- E) Build count matrix (features are "orig_seq + correct_annotation")
count_mat <- make_count_matrix(sample_tbls_with_genome)
saveRDS(count_mat, file.path(cfg$out_dir, "count_matrix.rds"))

# --- F) CPM + filtering
flt <- cpm_filter(count_mat, cfg$cpm_threshold, cfg$cpm_min_samples)
saveRDS(flt$cpm,  file.path(cfg$out_dir, "cpm_matrix.rds"))
saveRDS(flt$keep, file.path(cfg$out_dir, "cpm_keep_mask.rds"))

# Save filtered matrices (optional)
count_mat_keep <- count_mat[flt$keep, , drop = FALSE]
cpm_keep <- flt$cpm[flt$keep, , drop = FALSE]
saveRDS(count_mat_keep, file.path(cfg$out_dir, "count_matrix_cpm_filtered.rds"))
saveRDS(cpm_keep,       file.path(cfg$out_dir, "cpm_matrix_filtered.rds"))

# --- G) Minimal QC plots 
# 1) Composition per "set" (ELAC3S/GFP3S/WT3S/...)
plot_composition <- function(sample_tables, out_png) {
  dt <- data.table::rbindlist(sample_tables, use.names = TRUE, fill = TRUE)
  dt_sum <- dt[, .N, by = .(set, correct_RNA_type)]
  dt_sum[, perc := 100 * N / sum(N), by = set]
  
  p <- ggplot2::ggplot(
    dt_sum,
    ggplot2::aes(x = "", y = perc, fill = correct_RNA_type)
  ) +
    ggplot2::geom_col(color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~set, nrow = 2) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(out_png, p, width = 12, height = 6, dpi = 300)
}

plot_composition(
  sample_tbls_with_genome,
  file.path(cfg$out_dir, "rna_type_composition_by_set.png")
)

message("Done. Outputs written to: ", cfg$out_dir)