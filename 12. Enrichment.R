# ===== Enrichment analysis for Canis familiaris =====
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install needed packages
for (pkg in c("clusterProfiler", "org.Cf.eg.db", "AnnotationDbi", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# ---- SETTINGS ----
deg_dir <- "DEA2/annotated_DEGs"          # folder with annotated DEG CSVs
out_dir <- file.path(deg_dir, "enrichment_results")
dir.create(out_dir, showWarnings = FALSE)

# ---- Function to run enrichment for one file ----
run_enrichment <- function(file) {
  # Read file
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Use SYMBOL if available, otherwise PROBEID → map to SYMBOL
  if ("SYMBOL" %in% names(df)) {
    gene_symbols <- na.omit(unique(df$SYMBOL))
  } else {
    stop("No SYMBOL column in file: ", file)
  }
  
  # Map SYMBOL to Entrez IDs for enrichment
  gene_df <- bitr(gene_symbols, fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Cf.eg.db)
  
  entrez_ids <- unique(gene_df$ENTREZID)
  
  if (length(entrez_ids) < 5) {
    cat("⚠️ Skipping", basename(file), "— not enough genes for enrichment.\n")
    return(NULL)
  }
  
  # GO enrichment
  ego <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Cf.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "ALL",       # BP, CC, MF or ALL
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  
  # KEGG enrichment
  ekegg <- enrichKEGG(gene         = entrez_ids,
                      organism     = 'cfa',
                      pvalueCutoff = 0.05)
  
  # Save results
  base <- tools::file_path_sans_ext(basename(file))
  write.csv(as.data.frame(ego),
            file.path(out_dir, paste0(base, "_GO_enrichment.csv")),
            row.names = FALSE)
  write.csv(as.data.frame(ekegg),
            file.path(out_dir, paste0(base, "_KEGG_enrichment.csv")),
            row.names = FALSE)
  
  cat("✅ Enrichment saved for", basename(file), "\n")
}

# ---- Run for all DEG files ----
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE)
for (f in deg_files) {
  run_enrichment(f)
}

cat("\n🎯 All enrichment results saved in:", normalizePath(out_dir), "\n")


