# ===== Enrichment Visualization Suite — works with *_GO_enrichment.csv and *_KEGG_enrichment.csv =====
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Cf.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(AnnotationDbi)
  if (!requireNamespace("GOSemSim", quietly = TRUE)) BiocManager::install("GOSemSim", ask = FALSE)
  library(GOSemSim)
})

ENRICH_DIR <- "DEA2/annotated_DEGs/enrichment_results"
PLOT_DIR   <- file.path(ENRICH_DIR, "plots")
PLOT_GO    <- file.path(PLOT_DIR, "GO")
PLOT_KEGG  <- file.path(PLOT_DIR, "KEGG")
dir.create(PLOT_GO, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_KEGG, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
safe_gsave <- function(plot, path, w = 10, h = 7, dpi = 300) {
  try(ggsave(path, plot, width = w, height = h, dpi = dpi), silent = TRUE)
}

.semBP <- tryCatch(godata(annoDb = "org.Cf.eg.db", ont = "BP"), error = function(e) NULL)
has_edges <- function(er_or_sim) {
  ts <- tryCatch(er_or_sim@termsim, error = function(e) NULL)
  if (is.null(ts) || !is.matrix(ts) || nrow(ts) < 2) return(FALSE)
  any(ts[upper.tri(ts)] > 0, na.rm = TRUE)
}

load_as_enrichResult <- function(csv_file, type = c("GO","KEGG")) {
  df <- read.csv(csv_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(df)) return(NULL)
  type <- match.arg(type)
  if (!"p.adjust" %in% names(df)) {
    if ("qvalue" %in% names(df)) df$p.adjust <- df$qvalue
    else if ("pvalue" %in% names(df)) df$p.adjust <- p.adjust(df$pvalue, "BH")
    else df$p.adjust <- NA_real_
  }
  if (!"Description" %in% names(df)) {
    ch <- names(df)[sapply(df, is.character)]
    if (length(ch)) names(df)[names(df) == ch[1]] <- "Description"
  }
  if (!"geneID" %in% names(df)) df$geneID <- ""
  new("enrichResult",
      result        = df,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.2,
      organism      = if (type == "GO") "Canis familiaris" else "cfa",
      keytype       = "ENTREZID")
}

safe_emap_GO <- function(er, out_png, show_n = 30, padj_cut = 0.2) {
  if (is.null(er) || nrow(as.data.frame(er)) < 2 || is.null(.semBP)) return(FALSE)
  er2 <- tryCatch(simplify(er, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) er)
  er2 <- head(er2, n = show_n)
  er2 <- tryCatch(er2[er2@result$p.adjust < padj_cut, ], error = function(e) er2)
  if (nrow(as.data.frame(er2)) < 2) return(FALSE)
  sim <- try(pairwise_termsim(er2, semData = .semBP, method = "Wang"), silent = TRUE)
  if (inherits(sim, "try-error")) return(FALSE)
  if (!has_edges(sim)) return(FALSE)
  p <- try(emapplot(sim, showCategory = min(show_n, nrow(as.data.frame(sim))), layout = "kk"), silent = TRUE)
  if (inherits(p, "try-error")) return(FALSE)
  safe_gsave(p, out_png, w = 12, h = 9)
  TRUE
}

safe_emap_KEGG <- function(er, out_png, show_n = 30, padj_cut = 0.2) {
  if (is.null(er) || nrow(as.data.frame(er)) < 2) return(FALSE)
  er2 <- head(er, n = show_n)
  er2 <- tryCatch(er2[er2@result$p.adjust < padj_cut, ], error = function(e) er2)
  if (nrow(as.data.frame(er2)) < 2) return(FALSE)
  sim <- try(pairwise_termsim(er2), silent = TRUE)
  if (inherits(sim, "try-error")) return(FALSE)
  if (!has_edges(sim)) return(FALSE)
  p <- try(emapplot(sim, showCategory = min(show_n, nrow(as.data.frame(sim))), layout = "kk"), silent = TRUE)
  if (inherits(p, "try-error")) return(FALSE)
  safe_gsave(p, out_png, w = 12, h = 9)
  TRUE
}

make_plots <- function(res, tag, is_go = TRUE) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(invisible(NULL))
  topN <- 15
  out_base <- file.path(if (is_go) PLOT_GO else PLOT_KEGG, tag)
  
  p1 <- dotplot(res, showCategory = topN) + ggtitle(paste(tag, "— Dotplot"))
  safe_gsave(p1, paste0(out_base, "_dotplot.png"), w = 9, h = 7)
  
  p2 <- barplot(res, showCategory = topN) + ggtitle(paste(tag, "— Barplot"))
  safe_gsave(p2, paste0(out_base, "_barplot.png"), w = 9, h = 7)
  
  if (nrow(as.data.frame(res)) > 1) {
    p4 <- try(cnetplot(res, showCategory = min(topN, nrow(as.data.frame(res))), circular = FALSE, colorEdge = TRUE) +
                ggtitle(paste(tag, "— Cnetplot")), silent = TRUE)
    if (!inherits(p4, "try-error")) safe_gsave(p4, paste0(out_base, "_cnetplot.png"), w = 12, h = 9)
  }
  
  df <- as.data.frame(res)
  if (is_go && "ONTOLOGY" %in% names(df)) {
    p5 <- try(treeplot(res, showCategory = topN) + ggtitle(paste(tag, "— GO Treeplot")), silent = TRUE)
    if (!inherits(p5, "try-error")) safe_gsave(p5, paste0(out_base, "_treeplot.png"), w = 9, h = 7)
  }
  
  ok <- if (is_go) {
    safe_emap_GO(res, paste0(out_base, "_emapplot.png"), show_n = 30, padj_cut = 0.2)
  } else {
    safe_emap_KEGG(res, paste0(out_base, "_emapplot.png"), show_n = 30, padj_cut = 0.2)
  }
  if (!ok) {
    note <- if (is_go)
      "No GO term similarity edges (or GOSemSim unavailable). Dot/bar saved instead."
    else
      "No KEGG pathway overlaps among top terms. Dot/bar saved instead."
    writeLines(note, paste0(out_base, "_emapplot_NOTE.txt"))
  }
}

# ---------- main ----------
csv_files <- list.files(ENRICH_DIR,
                        pattern = "_(GO|KEGG)_enrichment\\.csv$",
                        full.names = TRUE)
cat("Found", length(csv_files), "enrichment CSVs in", ENRICH_DIR, "\n")
if (!length(csv_files)) stop("No enrichment CSVs found. Expected *_GO_enrichment.csv or *_KEGG_enrichment.csv")

for (f in csv_files) {
  tag  <- sub("\\.csv$", "", basename(f))
  is_go <- grepl("_GO_", f)
  type  <- if (is_go) "GO" else "KEGG"
  cat("Processing", tag, "...\n")
  er <- load_as_enrichResult(f, type)
  make_plots(er, tag, is_go = is_go)
}

cat("✅ Plots saved under:\n -", normalizePath(PLOT_GO), "\n -", normalizePath(PLOT_KEGG), "\n")