#########################################
# 1. Load packages
#########################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "limma", "canine2.db", "pheatmap", "ggplot2"
), ask = FALSE, update = FALSE)

library(limma)
library(canine2.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)

#########################################
# 2. File paths
#########################################
expr_file <- "/Users/alalehvd/Documents/GitHub/CAD-Analysis/expression_rma_probes_GSMcols.csv"
meta_file <- "/Users/alalehvd/Documents/GitHub/CAD-Analysis/samples_metadata.csv"
out_dir   <- "/Users/alalehvd/Documents/GitHub/CAD-Analysis"

#########################################
# 3. Load data
#########################################
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_file)

# Make sure columns match metadata order
expr <- expr[, meta$GSM]

#########################################
# 4. Map probes → gene symbols
#########################################
gene_symbols <- mapIds(
  canine2.db,
  keys = rownames(expr),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

expr_annot <- data.frame(GeneSymbol = gene_symbols, expr)
write.csv(expr_annot, file.path(out_dir, "expression_rma_genes.csv"), row.names = TRUE)

#########################################
# 5. Prepare design matrix (fixed)
#########################################
# Re-factor from scratch with explicit levels
meta$group <- factor(as.character(meta$group), levels = c("Healthy","NLAD","ALAD"))

# Build a fresh design from the data frame (not the old object)
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

#########################################
# 6. Sanity checks & alignment
#########################################
# Ensure all GSMs in metadata exist as columns in expr
if (!all(meta$GSM %in% colnames(expr))) {
  missing <- setdiff(meta$GSM, colnames(expr))
  stop(sprintf("These GSMs are in metadata but not in expression: %s",
               paste(missing, collapse = ", ")))
}
# Reorder columns to match metadata order exactly
expr <- expr[, meta$GSM]

# Quick checks
print(table(meta$group))          # expect 13 Healthy, 13 NLAD, 13 ALAD
print(colSums(design))            # expect 13 13 13
stopifnot(qr(design)$rank == ncol(design))  # full rank

#########################################
# 7. Differential expression (limma)
#########################################
fit  <- limma::lmFit(expr, design)
cont <- limma::makeContrasts(
  ALADvsHealthy = ALAD - Healthy,
  NLADvsHealthy = NLAD - Healthy,
  levels = design
)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

# Save full result tables
alad_res <- limma::topTable(fit2, coef = "ALADvsHealthy", number = Inf, sort.by = "P")
nlad_res <- limma::topTable(fit2, coef = "NLADvsHealthy", number = Inf, sort.by = "P")

write.csv(alad_res, file.path(out_dir, "DEG_ALAD_vs_Healthy.csv"))
write.csv(nlad_res, file.path(out_dir, "DEG_NLAD_vs_Healthy.csv"))

# Simple plots
plot_volcano <- function(res, title, outpng) {
  df <- data.frame(logFC = res$logFC, P = res$P.Value)
  df$logP <- -log10(df$P)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = logP)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    ggplot2::labs(title = title, x = "log2 fold-change", y = "-log10 p") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(filename = file.path(out_dir, outpng), plot = p, width = 6, height = 5, dpi = 150)
}
plot_volcano(alad_res, "ALAD vs Healthy", "volcano_ALADvsHealthy.png")
plot_volcano(nlad_res, "NLAD vs Healthy", "volcano_NLADvsHealthy.png")

topN <- 50
top_ids <- rownames(alad_res)[seq_len(min(topN, nrow(alad_res)))]
mat_top <- expr[top_ids, , drop = FALSE]
ann <- data.frame(Group = meta$group); rownames(ann) <- meta$GSM
pheatmap::pheatmap(mat_top, scale = "row",
                   annotation_col = ann,
                   filename = file.path(out_dir, "heatmap_top50_ALADvsHealthy.png"),
                   width = 8, height = 10)

cat("✅ Done. Results in:\n",
    "- DEG_ALAD_vs_Healthy.csv\n",
    "- DEG_NLAD_vs_Healthy.csv\n",
    "- volcano_ALADvsHealthy.png\n",
    "- volcano_NLADvsHealthy.png\n",
    "- heatmap_top50_ALADvsHealthy.png\n", sep = "")
