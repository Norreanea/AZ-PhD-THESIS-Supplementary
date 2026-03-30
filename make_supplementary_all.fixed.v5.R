#!/usr/bin/env Rscript
# make_supplementary_all.R
#
# One unified generator for Supplementary/all:
# - Static site in docs/ (GitHub Pages ready)
# - Pages: Home, UMAP, Figures, Tables, Scripts
# - Tables: Open preview (HTML per XLSX) + Download original
# - UMAP: plotly HTML widgets (general + detailed) + activity pages (72 hpa)
#
# Run from Supplementary/all:
#   Rscript make_supplementary_all.R
# or:
#   Rscript make_supplementary_all.R /path/to/Supplementary/all /path/to/AZ_THESIS_v2.6_PJ.docx

args <- commandArgs(trailingOnly = TRUE)
# ROOT can be provided as arg1. If not provided, fall back to a local default if it exists.
ROOT <- if (length(args) >= 1) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  default_root <- "G:/PhD_final/thesis/Supplementary_Interactive"
  if (dir.exists(default_root)) normalizePath(default_root, winslash = "/", mustWork = TRUE) else getwd()
}
# -------------------------
# CONFIG (edit if needed)
# -------------------------
DEFAULT_REDUCTION <- "umap.d33.nn100.md0.3"
COND_COL <- "condition_pretty"
POP_GENERAL_COL <- "general_population"
POP_DETAILED_COL <- "final_population"

# Inputs for UMAP static pages (recommended to place in ROOT/data/)
OBJ_PATH <- file.path(ROOT, "data", "result_obj.qs")          # Seurat qs
ACTIVITY_PATH <- file.path(ROOT, "data", "srna_activity.RData") # contains "best" (your activity table)

# Optional: thesis docx for extracting mentions of Supplementary Table S#
THESIS_DOCX <- if (length(args) >= 2) {
  normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
} else {
  default_docx <- "G:/PhD_final/thesis/AZaremba_THESIS-FINAL.docx"
  if (file.exists(default_docx)) {
    normalizePath(default_docx, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(ROOT, "..", "..", "AZaremba_THESIS-FINAL.docx"), winslash = "/", mustWork = FALSE)
  }
}
# UMAP generation controls
DOWNSAMPLE_N <- 0L                 # 0 = all cells; otherwise sample for smaller HTML
MAX_SIDS <- Inf                    # number of sRNAs to render (Inf = all)
EMBED_FIRST_N_ACTIVITY <- 10L      # embed first N activity iframes into UMAP index; rest as links

# Output
REPORT_DIR <- file.path(ROOT, "report")
DOCS_DIR   <- file.path(ROOT, "docs")
UMAP_STAGE_DIR <- file.path(REPORT_DIR, "umap_static")   # stage widgets outside docs/ (render_site may wipe docs/)
UMAP_DIR   <- UMAP_STAGE_DIR

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(UMAP_DIR,   recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Dependencies
# -------------------------
need <- c("rmarkdown", "DT", "readxl", "officer", "qs", "Seurat", "plotly", "htmlwidgets", "htmltools")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall: install.packages(c(", paste(sprintf('\"%s\"', missing), collapse = ", "), "))")
}

# -------------------------
# Small helpers
# -------------------------
write_lines <- function(path, x) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, con = path, useBytes = TRUE)
}

copy_dir_recursive <- function(from, to) {
  if (!dir.exists(from)) return(invisible(NULL))
  if (dir.exists(to)) unlink(to, recursive = TRUE, force = TRUE)
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(from, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (length(files) == 0) return(invisible(NULL))
  for (f in files) {
    rel <- substring(f, nchar(from) + 2)
    dest <- file.path(to, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, dest, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
  invisible(NULL)
}

slugify <- function(x) gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(x))
href_file <- function(dir, fname) paste0(dir, "/", utils::URLencode(fname, reserved = FALSE))

# Minimal string concatenation helper (keeps long JS blocks readable).
`%+%` <- function(a, b) paste0(a, b)

# Clean Word/Zotero field-codes from extracted thesis sentences.
clean_mention_text <- function(x) {
  if (length(x) == 0) return(character())
  x <- as.character(x)
  x <- gsub("\\r\\n|\\r|\\n", " ", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  # Remove Zotero CSL_CITATION field blobs that officer can expose as plain text.
  x <- gsub("ADDIN\\s+ZOTERO_ITEM\\s+CSL_CITATION\\s*\\{.*\\}", "", x, perl = TRUE)
  x <- gsub("ADDIN\\s+ZOTERO_ITEM\\s+CSL_CITATION", "", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

# -------------------------
# Palettes (paste your exact vectors)
# -------------------------
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

DETAILED_COLS <- c(
  "Early epidermal progenitor"             = "#ABCEDF",
  "Late epidermal progenitor"              = "#6693C4",
  "Epidermis (broad)"                      = "#5685BD",
  "Epidermis (multiciliated)"              = "#3468B0",
  "Epidermal progenitor (multiciliated)"   = "#4576B6",
  "Epidermal secretory gland progenitor"   = "#99BFD8",
  "Eye progenitor"                         = "#00D3D3",
  "Basal cell"                             = "#D3E4BA",
  "Goblet cell"                            = "#8DB180",
  "Phagocyte (broad)"                      = "#497F46",
  "Body pigment progenitor"                = "#C08A5A",
  "Body pigment cell"                      = "#8E5A3C",
  "Muscle progenitor"                      = "#F2C7CA",
  "BWM (dorsal midline)"                   = "#E6A3A6",
  "ECM-producing muscle"                   = "#CC5C5D",
  "Posterior pole/PCG muscle"              = "#C1393A",
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
  "Pharyngeal epithelium"                  = "#FFFF00",
  "Pharyngeal progenitor"                  = "#D9D900",
  "Pharyngeal phagocytic-type cell"        = "#C9C900",
  "Protonephridial flame cell"             = "#874A68",
  "Protonephridial tubule cell"            = "#87345F"
)

# -------------------------
# Thesis mention extraction (optional)
# -------------------------
extract_supp_table_mentions <- function(docx_path) {
  if (!file.exists(docx_path)) return(list())
  doc <- officer::read_docx(docx_path)
  s   <- officer::docx_summary(doc)
  txt <- s$text
  txt <- txt[nzchar(trimws(txt))]
  txt <- clean_mention_text(txt)
  out <- list()
  for (k in 1:60) {
    key <- paste0("S", k)
    pat <- sprintf("(Supplementary|Suplementary)\\s+Table\\s+S%s\\b", k)
    hits <- unique(grep(pat, txt, ignore.case = TRUE, value = TRUE))
    if (length(hits) > 0) out[[key]] <- hits
  }
  out
}
MENTIONS <- extract_supp_table_mentions(THESIS_DOCX)

# -------------------------
# Curated table descriptions (edit freely)
# Keys must match filenames in tables/
# -------------------------
TABLE_DESCRIPTIONS <- list(
  "Table S1.xlsx"  = "Dataset-level metadata for the sequencing datasets used in this thesis.",
  "Table S2.xlsx"  = "Bulk small RNA-seq support table used for sncRNA differential analysis and summaries.",
  "Table S3.xlsx"  = "Bulk long RNA-seq support table used for differential gene expression and downstream GO analyses.",
  "Table S4.xlsx"  = "Curated marker reference used for automated cell-type assignment and validation of population labels.",
  "Table S5.xlsx"  = "Final population labels used across the scRNA-seq atlas analyses.",
  "Table S6.xlsx"  = "Population-resolved scRNA-seq differential expression signatures and summary statistics.",
  "Table S7.xlsx"  = "Δ (delta) connectivity measures between Smed ELAC2 KD and controls per edge and time point.",
  "Table S8.xlsx"  = "sncRNA raw reads used for low-input sncRNA-seq downstream comparisons.",
  "Table S9.xlsx"  = "sncRNA activity inferred from scRNA-seq of miRNAs, piRNAs and tRFs that were significantly affected by Smed ELAC2 silencing in bulk RNA-seq."
)

# -------------------------
# 1) Generate static UMAP HTML widgets into docs/umap_static/
# -------------------------
generate_umap_static <- function() {
  if (!file.exists(OBJ_PATH)) {
    message("UMAP inputs not found: ", OBJ_PATH, " (UMAP section will still render but without widgets).")
    return(invisible(list(ok = FALSE, pop_general_file = NA_character_, pop_detailed_file = NA_character_, activity_map = NULL)))
  }
  
  obj <- qs::qread(OBJ_PATH)
  
  if (!DEFAULT_REDUCTION %in% names(obj@reductions)) {
    stop("Reduction not found: ", DEFAULT_REDUCTION)
  }
  need_meta <- c(COND_COL, POP_GENERAL_COL, POP_DETAILED_COL)
  if (!all(need_meta %in% colnames(obj@meta.data))) {
    stop("Missing meta columns in Seurat object: ", paste(setdiff(need_meta, colnames(obj@meta.data)), collapse = ", "))
  }
  
  emb <- Seurat::Embeddings(obj, reduction = DEFAULT_REDUCTION)[, 1:2, drop = FALSE]
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  
  df <- data.frame(
    cell = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    condition = as.character(obj@meta.data[rownames(emb), COND_COL]),
    pop_general  = as.character(obj@meta.data[rownames(emb), POP_GENERAL_COL]),
    pop_detailed = as.character(obj@meta.data[rownames(emb), POP_DETAILED_COL]),
    stringsAsFactors = FALSE
  )
  
  if (is.numeric(DOWNSAMPLE_N) && DOWNSAMPLE_N > 0 && nrow(df) > DOWNSAMPLE_N) {
    set.seed(1)
    df <- df[sample.int(nrow(df), DOWNSAMPLE_N), , drop = FALSE]
  }
  
  all_conditions <- sort(unique(df$condition))
  cond72 <- all_conditions[grepl("72hpa", all_conditions, ignore.case = TRUE)]
  cond72 <- unique(c(cond72[grepl("^WT", cond72)], cond72[grepl("^GFP", cond72)], cond72[grepl("^ELAC", cond72)]))
  
  plot_umap_pop <- function(df_use, pop_col, pal, title) {
    pop <- as.character(df_use[[pop_col]])
    col_hex <- pal[pop]
    col_hex[is.na(col_hex)] <- "#BDBDBD"
    hover <- paste0(
      "Cell: ", df_use$cell,
      "<br>Condition: ", df_use$condition,
      "<br>General: ", df_use$pop_general,
      "<br>Detailed: ", df_use$pop_detailed
    )
    plotly::plot_ly(
      data = df_use,
      x = ~UMAP_1, y = ~UMAP_2,
      type = "scattergl", mode = "markers",
      text = hover, hoverinfo = "text",
      marker = list(size = 1.6, color = col_hex, opacity = 0.9),
      showlegend = FALSE
    ) |>
      plotly::layout(
        title = list(text = title, x = 0.02, xanchor = "left"),
        xaxis = list(title = "", zeroline = FALSE),
        yaxis = list(title = "", zeroline = FALSE),
        margin = list(l = 10, r = 10, b = 10, t = 40)
      ) |>
      plotly::config(responsive = TRUE, displaylogo = FALSE)
  }
  
  # Save (shared libs)
  save_widget <- function(w, filename) {
    htmlwidgets::saveWidget(
      widget = w,
      file = file.path(UMAP_DIR, filename),
      selfcontained = FALSE,
      # IMPORTANT: keep libdir RELATIVE so widgets work after moving docs/ or on GitHub Pages.
      libdir = "libs"
    )
    filename
  }
  
  pop_general_file  <- save_widget(plot_umap_pop(df, "pop_general",  FAMILY_COLORS,  "UMAP – general populations"),
                                   "pop_general.html")
  pop_detailed_file <- save_widget(plot_umap_pop(df, "pop_detailed", DETAILED_COLS, "UMAP – detailed populations"),
                                   "pop_detailed.html")
  
  # Activity handling (expects your "best" long table)
  load_activity_object <- function(rdata_path) {
    if (!file.exists(rdata_path)) return(NULL)
    env <- new.env(parent = emptyenv())
    base::load(rdata_path, envir = env)
    objs <- ls(env, all.names = TRUE)
    if (length(objs) == 0) return(NULL)
    preferred <- c("best", "activity", "activity_obj", "srna_activity")
    pick <- preferred[preferred %in% objs]
    name <- if (length(pick) >= 1) pick[[1]] else objs[[1]]
    env[[name]]
  }
  
  activity_obj <- load_activity_object(ACTIVITY_PATH)
  activity_map <- data.frame(sid = character(), file = character(), stringsAsFactors = FALSE)
  
  build_activity_from_best <- function(best_df, df72, sid) {
    best_df <- as.data.frame(best_df)
    stopifnot(all(c("sRNA", "celltype_use", "WT", "GFP", "ELAC") %in% colnames(best_df)))
    
    sub <- best_df[as.character(best_df$sRNA) == sid, , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    
    if ("timepoint" %in% colnames(sub)) {
      sub <- sub[as.integer(sub$timepoint) == 72, , drop = FALSE]
      if (nrow(sub) == 0) return(NULL)
    }
    
    if ("model" %in% colnames(sub) && length(unique(sub$model)) > 1) {
      if ("score_expected" %in% colnames(sub) && is.numeric(sub$score_expected)) {
        model_pick <- sub$model[which.max(sub$score_expected)][1]
      } else {
        model_pick <- sub$model[1]
      }
      sub <- sub[as.character(sub$model) == as.character(model_pick), , drop = FALSE]
    }
    
    lookup <- sub[, c("celltype_use", "WT", "GFP", "ELAC"), drop = FALSE]
    
    cond_key <- ifelse(grepl("^WT", df72$condition), "WT",
                       ifelse(grepl("^GFP", df72$condition), "GFP",
                              ifelse(grepl("^ELAC", df72$condition), "ELAC", NA_character_)))
    
    idx <- match(df72$pop_detailed, lookup$celltype_use)
    
    out <- rep(NA_real_, nrow(df72))
    out[cond_key == "WT"]   <- as.numeric(lookup$WT[idx[cond_key == "WT"]])
    out[cond_key == "GFP"]  <- as.numeric(lookup$GFP[idx[cond_key == "GFP"]])
    out[cond_key == "ELAC"] <- as.numeric(lookup$ELAC[idx[cond_key == "ELAC"]])
    names(out) <- df72$cell
    out
  }
  
  plot_activity_72 <- function(df72, act_vec, sid, cond72) {
    d <- df72[df72$condition %in% cond72, , drop = FALSE]
    d$activity <- as.numeric(act_vec[d$cell])
    
    vals <- d$activity[is.finite(d$activity)]
    lim <- if (length(vals) == 0) 1 else max(abs(stats::quantile(vals, probs = c(0.10, 0.90), na.rm = TRUE)))
    if (!is.finite(lim) || lim == 0) lim <- 1
    
    cs <- list(list(0, "#0000FF"), list(0.5, "#FFFFFF"), list(1, "#FF0000"))
    
    make_one <- function(cond_label, show_scale = FALSE) {
      dd <- d[d$condition == cond_label, , drop = FALSE]
      hover <- paste0(
        "Cell: ", dd$cell,
        "<br>Condition: ", dd$condition,
        "<br>Detailed: ", dd$pop_detailed,
        "<br>", sid,
        "<br>Activity: ", signif(dd$activity, 4)
      )
      plotly::plot_ly(
        data = dd,
        x = ~UMAP_1, y = ~UMAP_2,
        type = "scattergl", mode = "markers",
        text = hover, hoverinfo = "text",
        marker = list(
          size = 1.6,
          color = ~activity,
          colorscale = cs,
          cmin = -lim, cmax = lim,
          showscale = isTRUE(show_scale),
          colorbar = list(title = "Activity", len = 0.85),
          opacity = 0.9
        ),
        showlegend = FALSE
      ) |>
        plotly::layout(
          title = list(text = gsub("([0-9]+hpa)", "<br>\\1", cond_label, perl = TRUE),
                       x = 0.02, xanchor = "left"),
          xaxis = list(title = "", zeroline = FALSE),
          yaxis = list(title = "", zeroline = FALSE),
          margin = list(l = 10, r = 10, b = 10, t = 40)
        ) |>
        plotly::config(responsive = TRUE, displaylogo = FALSE)
    }
    
    ord <- unique(c(cond72[grepl("^WT", cond72)], cond72[grepl("^GFP", cond72)], cond72[grepl("^ELAC", cond72)]))
    ord <- ord[ord %in% unique(d$condition)]
    if (length(ord) == 0) ord <- unique(d$condition)
    
    plist <- lapply(seq_along(ord), function(i) make_one(ord[[i]], show_scale = (i == length(ord))))
    
    plotly::subplot(plist, nrows = 1, shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE) |>
      plotly::layout(
        title = list(text = sid, x = 0.02, xanchor = "left"),
        margin = list(l = 10, r = 10, b = 10, t = 50)
      ) |>
      plotly::config(responsive = TRUE, displaylogo = FALSE)
  }
  
  if (!is.null(activity_obj) && is.data.frame(activity_obj) && length(cond72) > 0) {
    df72 <- df[df$condition %in% cond72, , drop = FALSE]
    sids <- sort(unique(as.character(activity_obj$sRNA)))
    if (is.finite(MAX_SIDS)) sids <- head(sids, MAX_SIDS)
    
    for (i in seq_along(sids)) {
      sid <- sids[[i]]
      act_vec <- build_activity_from_best(activity_obj, df72 = df72, sid = sid)
      if (is.null(act_vec)) next
      
      w <- plot_activity_72(df72, act_vec, sid = sid, cond72 = cond72)
      fname <- paste0(slugify(paste0("activity__", sid)), ".html")
      fname <- save_widget(w, fname)
      activity_map <- rbind(activity_map, data.frame(sid = sid, file = fname, stringsAsFactors = FALSE))
    }
  }
  
  # UMAP index.html (static)
  css <- "
  body { font-family: -apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif; margin: 24px; }
  iframe { width: 100%; height: 560px; border: 1px solid #ddd; border-radius: 6px; }
  .grid2 { display:grid; grid-template-columns: 1fr 1fr; gap: 14px; }
  .small { color:#555; font-size:0.95rem; }
  .code { font-family: ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace; }
  .controls { display:flex; gap: 10px; align-items:center; flex-wrap: wrap; margin: 8px 0 10px; }
  .controls input { min-width: 520px; max-width: 100%; padding: 6px 8px; }
  "

  # Write activity index for other pages (and for debugging).
  if (!is.null(activity_map) && nrow(activity_map) > 0) {
    utils::write.table(activity_map, file = file.path(UMAP_DIR, "activity_index.tsv"),
                       sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
    utils::write.table(activity_map, file = file.path(REPORT_DIR, "activity_index.tsv"),
                       sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
  }

  activity_block <- if (is.null(activity_map) || nrow(activity_map) == 0) {
    htmltools::tags$p(class = "small", "No activity pages generated (missing activity input or no 72 hpa conditions).")
  } else {
    # Use datalist for fast searching even with many sRNAs.
    first_sid  <- activity_map$sid[[1]]
    first_file <- activity_map$file[[1]]
    htmltools::tagList(
      htmltools::tags$p(class = "small",
                        "Type/select an sRNA to load its 72 hpa activity UMAP (interactive; pan/zoom/hover)."),
      htmltools::tags$div(class = "controls",
        htmltools::tags$label("sRNA:", `for` = "activity_input"),
        htmltools::tags$input(id = "activity_input", list = "activity_list", value = first_sid),
        htmltools::tags$a(id = "activity_open", href = first_file, target = "_blank", rel = "noopener", "Open in new tab")
      ),
      htmltools::tags$datalist(
        id = "activity_list",
        lapply(seq_len(nrow(activity_map)), function(i) {
          htmltools::tags$option(value = activity_map$sid[[i]], `data-file` = activity_map$file[[i]])
        })
      ),
      htmltools::tags$iframe(id = "activity_frame", src = first_file, loading = "lazy"),
      htmltools::tags$script(htmltools::HTML(
        "(function(){\n" %+%
          "function fileForSid(sid){\n" %+%
            "var opts=document.querySelectorAll('#activity_list option');\n" %+%
            "for(var i=0;i<opts.length;i++){ if(opts[i].value===sid){ return opts[i].getAttribute('data-file'); } }\n" %+%
            "return null;\n" %+%
          "}\n" %+%
          "var input=document.getElementById('activity_input');\n" %+%
          "// Convenience: allow changing selection without manually deleting the current value\n" %+%
          "input.addEventListener('focus', function(){ try{ this.select(); }catch(e){} });\n" %+%
          "input.addEventListener('click', function(){ try{ this.select(); }catch(e){} });\n" %+%
          "var frame=document.getElementById('activity_frame');\n" %+%
          "var open=document.getElementById('activity_open');\n" %+%
          "function update(){\n" %+%
            "var f=fileForSid(input.value);\n" %+%
            "if(f){ frame.src=f; open.href=f; }\n" %+%
          "}\n" %+%
          "input.addEventListener('change', update);\n" %+%
          "input.addEventListener('input', update);\n" %+%
        "})();\n"
      ))
    )
  }
  index <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title("UMAP – interactive"),
      htmltools::tags$style(htmltools::HTML(css))
    ),
    htmltools::tags$body(
      htmltools::tags$h1("UMAP – interactive"),
      htmltools::tags$p(class = "small",
                        "UMAP embedding: ", htmltools::tags$span(class = "code", DEFAULT_REDUCTION), ". ",
                        "Population UMAPs are fully interactive (pan/zoom/hover). Activity pages are per sRNA (72 hpa only)."
      ),
      htmltools::tags$h2("Population UMAPs"),
      htmltools::tags$div(class = "grid2",
                          htmltools::tags$iframe(src = pop_general_file, loading = "lazy"),
                          htmltools::tags$iframe(src = pop_detailed_file, loading = "lazy")
      ),
      htmltools::tags$h2("Small RNA activity (72 hpa)"),
      activity_block
    )
  )
  
  htmltools::save_html(index, file = file.path(UMAP_DIR, "index.html"))
  invisible(list(ok = TRUE, pop_general_file = pop_general_file, pop_detailed_file = pop_detailed_file, activity_map = activity_map))
}

UMAP_INFO <- generate_umap_static()

# -------------------------
# 2) Create the R Markdown website pages (report/)
# -------------------------
site_yml <- c(
  "name: \"Supplementary\"",
  "output_dir: \"../docs\"",
  "navbar:",
  "  title: \"Supplementary\"",
  "  left:",
  "    - text: \"Home\"",
  "      href: index.html",
  "    - text: \"UMAP\"",
  "      href: umap.html",
  "    - text: \"Figures\"",
  "      href: figures.html",
  "    - text: \"Tables\"",
  "      href: tables.html",
  "    - text: \"Scripts\"",
  "      href: scripts.html",
  "output:",
  "  html_document:",
  "    theme: flatly",
  "    toc: true",
  "    toc_float: true",
  "    df_print: paged"
)
write_lines(file.path(REPORT_DIR, "_site.yml"), site_yml)

index_rmd <- c(
  "---",
  "title: \"Supplementary material – interactive report\"",
  "output: html_document",
  "---",
  "",
  "This page contains the supplementary materials of the thesis:",
  "",
  "- **UMAP** – interactive population maps and per-sRNA activity pages (72 hpa only).",
  "- **Figures** – PDFs (Open/Download).",
  "- **Tables** – XLSX tables rendered as interactive previews; original XLSX downloadable.",
  "- **Scripts** – code files (Open/Download) and previews for `.R` and `.sh`."
  
  
)
write_lines(file.path(REPORT_DIR, "index.Rmd"), index_rmd)

# ---- UMAP page (embeds static UMAP index via iframe) ----
# ---- UMAP page (embeds static UMAP index via iframe) ----
umap_rmd <- c(
  "---",
  "title: \"UMAP\"",
  "output: html_document",
  "---",
  "",
  "Population UMAPs are shown side-by-side. Small RNA activity UMAPs (72 hpa) can be loaded via a searchable selector.",
  "",
  if (is.list(UMAP_INFO) && isTRUE(UMAP_INFO$ok)) {
    c(
      "<style>",
      ".umap-grid{display:grid;grid-template-columns:1fr 1fr;gap:14px;}",
      ".umap-grid iframe{width:100%;height:560px;border:1px solid #ddd;border-radius:6px;}",
      ".act-controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin:8px 0 10px;}",
      ".act-controls input{min-width:520px;max-width:100%;padding:6px 8px;}",
      "</style>",
      "",
      "## Population UMAPs",
      "",
      "<div class=\"umap-grid\">",
      sprintf("<iframe src=\"umap_static/%s\" loading=\"lazy\"></iframe>", UMAP_INFO$pop_general_file),
      sprintf("<iframe src=\"umap_static/%s\" loading=\"lazy\"></iframe>", UMAP_INFO$pop_detailed_file),
      "</div>",
      "",
      "## Small RNA activity (72 hpa)",
      "",
      "```{r, echo=FALSE, message=FALSE, warning=FALSE}",
      "act_idx_path <- 'activity_index.tsv'",
      "if (file.exists(act_idx_path)) {",
      "  act <- utils::read.delim(act_idx_path, sep='\\t', stringsAsFactors=FALSE)",
      "  if (nrow(act) > 0) {",
      "    first_sid  <- act$sid[[1]]",
      "    first_file <- act$file[[1]]",
      "    js <- paste(",
      "      \"(function(){\",",
      "      \"  function fileForSid(sid){\",",
      "      \"    var opts=document.querySelectorAll('#activity_list option');\",",
      "      \"    for(var i=0;i<opts.length;i++){\",",
      "      \"      if(opts[i].value===sid){ return opts[i].getAttribute('data-file'); }\",",
      "      \"    }\",",
      "      \"    return null;\",",
      "      \"  }\",",
      "      \"  var input=document.getElementById('activity_input');\",",
      "      \"  var frame=document.getElementById('activity_frame');\",",
      "      \"  var open=document.getElementById('activity_open');\",",
      "      \"  if(!input.dataset.prev){ input.dataset.prev = input.value || ''; }\",",
      "      \"  function update(){\",",
      "      \"    var f=fileForSid(input.value);\",",
      "      \"    if(f){\",",
      "      \"      input.dataset.prev = input.value;\",",
      "      \"      frame.src='umap_static/'+f;\",",
      "      \"      open.href='umap_static/'+f;\",",
      "      \"    }\",",
      "      \"  }\",",
      "      \"  input.addEventListener('pointerdown', function(){\",",
      "      \"    if(!this.value) return;\",",
      "      \"    this.dataset.prev = this.value;\",",
      "      \"    this.value = '';\",",
      "      \"    this.setAttribute('placeholder', this.dataset.prev);\",",
      "      \"    setTimeout(function(){ input.dispatchEvent(new Event('input', {bubbles:true})); }, 0);\",",
      "      \"  });\",",
      "      \"  input.addEventListener('blur', function(){\",",
      "      \"    if(!this.value && this.dataset.prev){ this.value = this.dataset.prev; }\",",
      "      \"    this.removeAttribute('placeholder');\",",
      "      \"  });\",",
      "      \"  input.addEventListener('change', update);\",",
      "      \"  input.addEventListener('input', update);\",",
      "      \"})();\",",
      "      sep='\\n'",
      "    )",
      "    htmltools::tagList(",
      "      htmltools::tags$div(class='act-controls',",
      "        htmltools::tags$label('sRNA:', `for`='activity_input'),",
      "        htmltools::tags$input(",
      "          id='activity_input', type='text', list='activity_list', value=first_sid,",
      "          autocomplete='off', spellcheck='false'",
      "        ),",
      "        htmltools::tags$a(",
      "          id='activity_open', href=paste0('umap_static/', first_file),",
      "          target='_blank', rel='noopener', 'Open in new tab'",
      "        )",
      "      ),",
      "      htmltools::tags$datalist(id='activity_list',",
      "        lapply(seq_len(nrow(act)), function(i) {",
      "          htmltools::tags$option(value=act$sid[[i]], `data-file`=act$file[[i]])",
      "        })",
      "      ),",
      "      htmltools::tags$iframe(",
      "        id='activity_frame', src=paste0('umap_static/', first_file), loading='lazy',",
      "        style='width:100%;height:560px;border:1px solid #ddd;border-radius:6px;'",
      "      ),",
      "      htmltools::tags$script(htmltools::HTML(js))",
      "    )",
      "  } else {",
      "    cat('No activity pages found.')",
      "  }",
      "} else {",
      "  cat('Activity index not found (no activity pages generated).')",
      "}",
      "```",
      "",
      "(Full static page: [umap_static/index.html](umap_static/index.html))"
    )
  } else {
    "UMAP inputs were not found during generation."
  }
)
write_lines(file.path(REPORT_DIR, "umap.Rmd"), umap_rmd)

# ---- Figures page ----
figures_rmd <- c(
  "---",
  "title: \"Figures\"",
  "output: html_document",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
  "fig_dir <- normalizePath(file.path('..','figures'), winslash='/', mustWork=TRUE)",
  "fig_files <- sort(list.files(fig_dir, pattern='\\\\.pdf$', full.names=FALSE))",
  "```",
  "",
  "```{r}",
  "if (length(fig_files) == 0) {",
  "  cat('No PDF files found in figures/.\\n')",
  "} else {",
  "  href <- paste0('figures/', vapply(fig_files, utils::URLencode, character(1), reserved = FALSE))",
  "  df <- data.frame(",
  "    File = fig_files,",
  "    Open = sprintf('<a href=\"%s\" target=\"_blank\" rel=\"noopener\">Open</a>', href),",
  "    Download = sprintf('<a href=\"%s\" download>Download</a>', href),",
  "    stringsAsFactors = FALSE",
  "  )",
  "  DT::datatable(df, escape=FALSE, rownames=FALSE,",
  "    options=list(pageLength=25, scrollX=TRUE, dom='Bfrtip', buttons=c('copy','csv'))",
  "  )",
  "}",
  "```"
)
write_lines(file.path(REPORT_DIR, "figures.Rmd"), figures_rmd)

# ---- Scripts page ----
scripts_rmd <- c(
  "---",
  "title: \"Scripts\"",
  "output: html_document",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
  "scr_dir <- normalizePath(file.path('..','scripts'), winslash='/', mustWork=TRUE)",
  "scr_files <- sort(list.files(scr_dir, full.names=FALSE))",
  "```",
  "",
  "```{r}",
  "if (length(scr_files) == 0) {",
  "  cat('No files found in scripts/.\\n')",
  "} else {",
  "  ext <- tolower(tools::file_ext(scr_files))",
  "  href <- paste0('scripts/', vapply(scr_files, utils::URLencode, character(1), reserved = FALSE))",
  "  df <- data.frame(",
  "    File = scr_files,",
  "    Type = ifelse(ext == '', '(no extension)', ext),",
  "    Open = sprintf('<a href=\"%s\" target=\"_blank\" rel=\"noopener\">Open</a>', href),",
  "    Download = sprintf('<a href=\"%s\" download>Download</a>', href),",
  "    stringsAsFactors = FALSE",
  "  )",
  "  DT::datatable(df, escape=FALSE, rownames=FALSE,",
  "    options=list(pageLength=25, scrollX=TRUE, dom='Bfrtip', buttons=c('copy','csv'))",
  "  )",
  "}",
  "```",
  "",
  "## Previews (R and shell scripts)",
  "",
  "```{r, results='asis'}",
  "preview_one <- function(fname) {",
  "  path <- file.path(scr_dir, fname)",
  "  ext <- tolower(tools::file_ext(fname))",
  "  if (!file.exists(path)) return(invisible(NULL))",
  "  if (!ext %in% c('r','sh','bash')) return(invisible(NULL))",
  "  lang <- if (ext == 'r') 'r' else 'bash'",
  "  href <- paste0('scripts/', utils::URLencode(fname, reserved = FALSE))",
  "  cat('### ', fname, '\\n\\n', sep='')",
  "  cat(sprintf('- <a href=\"%s\" target=\"_blank\" rel=\"noopener\">Open</a>  \\\\n- <a href=\"%s\" download>Download</a>\\n\\n', href, href))",
  "  cat('```', lang, '\\n', sep='')",
  "  cat(paste(readLines(path, warn=FALSE), collapse='\\n'))",
  "  cat('\\n```\\n\\n')",
  "}",
  "for (f in scr_files) preview_one(f)",
  "```"
)
write_lines(file.path(REPORT_DIR, "scripts.Rmd"), scripts_rmd)

# -------------------------
# Tables: generate per-XLSX preview pages + main tables index
# -------------------------
TAB_SRC_DIR <- file.path(ROOT, "tables")
TAB_FILES <- sort(list.files(TAB_SRC_DIR, pattern = "\\.xlsx$", full.names = FALSE))

# Map a table filename to the expected Supplementary Table key (S1, S2, ...).
infer_s_key <- function(fn) {
  base <- sub("\\.xlsx$", "", fn)
  base <- sub("^Table[ _]?", "", base)
  base <- gsub("[^A-Za-z0-9]", "", base)
  if (grepl("^S[0-9]+$", base)) return(base)
  if (grepl("^[0-9]+$", base)) return(paste0("S", base))
  base
}

# Precompute cleaned thesis mentions for the Tables index page.
if (length(TAB_FILES) > 0) {
  tab_mentions <- data.frame(File = TAB_FILES, s_key = vapply(TAB_FILES, infer_s_key, character(1)),
                             Mention = "", stringsAsFactors = FALSE)
  tab_mentions$Mention <- vapply(tab_mentions$s_key, function(k) {
    q <- MENTIONS[[k]]
    if (is.null(q) || length(q) == 0) return("")
    q <- clean_mention_text(unique(q))
    q <- q[nzchar(q)]
    if (length(q) == 0) return("")
    paste(head(q, 2), collapse = " ")
  }, character(1))
  utils::write.table(tab_mentions[, c("File", "Mention"), drop = FALSE],
                     file = file.path(REPORT_DIR, "table_mentions.tsv"),
                     sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE)
}

# Per-table preview pages: table_<slug>.Rmd (renders to docs/table_<slug>.html)
for (fname in TAB_FILES) {
  slug <- slugify(fname)
  rmd_path <- file.path(REPORT_DIR, paste0("table_", slug, ".Rmd"))
  
  s_key <- infer_s_key(fname)
  
  desc <- TABLE_DESCRIPTIONS[[fname]]
  if (is.null(desc) || !nzchar(desc)) desc <- ""
  
  quotes <- MENTIONS[[s_key]]
  quote_block <- ""
  if (!is.null(quotes) && length(quotes) > 0) {
    quotes <- unique(quotes)
    quotes <- quotes[seq_len(min(length(quotes), 2))]
    quotes <- clean_mention_text(quotes)
    quotes <- quotes[nzchar(quotes)]
    if (length(quotes) > 0) {
      quote_block <- paste0(
        "\n\n**Mention in thesis (for traceability):**\n\n",
        paste(paste0("> ", quotes), collapse = "\n>\n"),
        "\n"
      )
    }
  }
  
  down_href <- href_file("tables", fname)
  
  rmd <- c(
    "---",
    paste0("title: \"", fname, "\""),
    "output: html_document",
    "---",
    "",
    if (nzchar(desc)) desc else NULL,
    if (nzchar(quote_block)) quote_block else NULL,
    "",
    paste0("- <a href=\"", down_href, "\" download>Download original XLSX</a>"),
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
    "tab_dir <- normalizePath(file.path('..','tables'), winslash='/', mustWork=TRUE)",
    paste0("fname <- ", deparse(fname)),
    "path_abs <- file.path(tab_dir, fname)",
    "```",
    "",
    "```{r}",
    "sheets <- tryCatch(readxl::excel_sheets(path_abs), error = function(e) character())",
    "if (length(sheets) == 0) {",
    "  htmltools::tags$em('(No readable sheets.)')",
    "} else {",
    "  blocks <- lapply(sheets, function(sh) {",
    "    dat <- tryCatch(readxl::read_excel(path_abs, sheet = sh), error = function(e) NULL)",
    "    if (is.null(dat)) {",
    "      return(htmltools::tagList(",
    "        htmltools::tags$h2(paste0('Sheet: ', sh)),",
    "        htmltools::tags$em('(Could not read this sheet.)')",
    "      ))",
    "    }",
    "    dt <- DT::datatable(as.data.frame(dat), rownames=FALSE, extensions=c('Buttons','Scroller'),",
    "      options=list(pageLength=25, scrollX=TRUE, scrollY=520, scroller=TRUE, dom='Bfrtip',",
    "        buttons=c('copy','csv','excel'), searchHighlight=TRUE))",
    "    htmltools::tagList(",
    "      htmltools::tags$h2(paste0('Sheet: ', sh)),",
    "      dt",
    "    )",
    "  })",
    "  htmltools::tagList(blocks)",
    "}",
    "```"
  )
  rmd <- rmd[!vapply(rmd, is.null, logical(1))]
  write_lines(rmd_path, rmd)
}

# Main tables index page (Open preview -> table_<slug>.html)
tables_rmd <- c(
  "---",
  "title: \"Tables\"",
  "output: html_document",
  "---",
  "",
  "Click **Open preview** to view an interactive HTML rendering of the XLSX file in the browser.",
  "Use **Download** to save the original `.xlsx` file.",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
  "tab_dir <- normalizePath(file.path('..','tables'), winslash='/', mustWork=TRUE)",
  "tab_files <- sort(list.files(tab_dir, pattern='\\\\.xlsx$', full.names=FALSE))",
  "mention_path <- file.path('table_mentions.tsv')",
  "mentions <- if (file.exists(mention_path)) utils::read.delim(mention_path, sep='\\t', stringsAsFactors=FALSE) else data.frame(File=character(), Mention=character())",
  "```",
  "",
  "```{r}",
  "if (length(tab_files) == 0) {",
  "  cat('No XLSX files found in tables/.\\n')",
  "} else {",
  "  slug <- gsub('[^A-Za-z0-9]+','_', tools::file_path_sans_ext(tab_files))",
  "  prev <- paste0('table_', slug, '.html')",
  "  down <- paste0('tables/', vapply(tab_files, utils::URLencode, character(1), reserved = FALSE))",
  "  m <- mentions$Mention[match(tab_files, mentions$File)]",
  "  m[is.na(m)] <- ''",
  "  m_html <- ifelse(nzchar(m), sprintf('<details><summary>Show</summary>%s</details>', htmltools::htmlEscape(m)), '')",
  "  df <- data.frame(",
  "    File = tab_files,",
  "    Mention = m_html,",
  "    Open = sprintf('<a href=\"%s\">Open preview</a>', prev),",
  "    Download = sprintf('<a href=\"%s\" download>Download</a>', down),",
  "    stringsAsFactors = FALSE",
  "  )",
  "  DT::datatable(df, escape=FALSE, rownames=FALSE,",
  "    options=list(pageLength=25, scrollX=TRUE, dom='Bfrtip', buttons=c('copy','csv'))",
  "  )",
  "}",
  "```"
)
write_lines(file.path(REPORT_DIR, "tables.Rmd"), tables_rmd)

# -------------------------
# 3) Render site -> docs/
# -------------------------
oldwd <- getwd()
setwd(REPORT_DIR)
on.exit(setwd(oldwd), add = TRUE)
rmarkdown::render_site(encoding = "UTF-8")

# -------------------------
# 4) Copy assets into docs/ so Open/Download work
# -------------------------
copy_dir_recursive(file.path(ROOT, "figures"), file.path(DOCS_DIR, "figures"))
copy_dir_recursive(file.path(ROOT, "scripts"), file.path(DOCS_DIR, "scripts"))
copy_dir_recursive(file.path(ROOT, "tables"),  file.path(DOCS_DIR, "tables"))
# Copy staged UMAP widgets into docs/ AFTER render_site() (render_site may wipe docs/)
copy_dir_recursive(UMAP_DIR, file.path(DOCS_DIR, "umap_static"))

message("Done.")
message("Open locally: ", normalizePath(file.path(DOCS_DIR, "index.html"), winslash = "/", mustWork = TRUE))

