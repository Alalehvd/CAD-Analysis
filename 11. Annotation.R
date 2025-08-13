# ===== Annotate all DEG CSV files =====
library(AnnotationDbi)

# ---- SETTINGS ----
deg_dir       <- "DEA2"          # folder with your DEG CSVs
out_dir       <- file.path(deg_dir, "annotated_DEGs")
platform_pkg  <- "canine2.db"    # replace if annotation(abatch) shows a different one

# ---- Setup ----
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!require(platform_pkg, character.only = TRUE)) {
  BiocManager::install(platform_pkg)
  library(platform_pkg, character.only = TRUE)
}
dir.create(out_dir, showWarnings = FALSE)

# ---- Function to annotate one file ----
annotate_deg <- function(deg_file, out_file) {
  deg <- read.csv(deg_file, row.names = 1, check.names = FALSE)
  probes <- rownames(deg)
  
  ann <- AnnotationDbi::select(get(platform_pkg),
                               keys = probes,
                               columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                               keytype = "PROBEID")
  
  # Collapse duplicates
  ann_collapsed <- aggregate(. ~ PROBEID, ann, function(x) paste(unique(na.omit(x)), collapse = "; "))
  
  # Merge
  deg$PROBEID <- rownames(deg)
  deg_annot <- merge(deg, ann_collapsed, by = "PROBEID", all.x = TRUE)
  
  write.csv(deg_annot, out_file, row.names = FALSE)
  cat("✅ Annotated file saved:", out_file, "\n")
}

# ---- Annotate all CSVs in the folder ----
deg_files <- list.files(deg_dir, pattern = "\\.csv$", full.names = TRUE)
for (f in deg_files) {
  out_path <- file.path(out_dir, basename(f))
  annotate_deg(f, out_path)
}

cat("\n🎯 All annotated files saved in:", normalizePath(out_dir), "\n")