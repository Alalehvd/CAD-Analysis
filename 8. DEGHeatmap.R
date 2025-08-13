# Heatmap (top N significant genes by FDR)
library(pheatmap)

expr <- read.csv("normalized_dataset.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("samples_metadata_grouped.csv", stringsAsFactors = FALSE)

# ensure GSM columns:
colnames(expr) <- sub("^(GSM\\d+).*", "\\1", colnames(expr))
meta <- meta[match(colnames(expr), meta$GSM), ]
ann <- data.frame(group = factor(meta$group, levels=c("Healthy","NLAD","ALAD")))
rownames(ann) <- colnames(expr)

heatmap_from_deg <- function(sig_file, tag, topN=50) {
  deg <- read.csv(sig_file, row.names = 1, check.names = FALSE)
  if (nrow(deg) == 0) return(invisible(NULL))
  sel <- head(rownames(deg[order(deg$adj.P.Val), , drop=FALSE]), topN)
  X <- as.matrix(expr[sel, , drop=FALSE])
  # z-score by gene
  Xz <- t(scale(t(X)))
  png(file.path("DEA/plots", paste0("heatmap_", tag, "_top", topN, ".png")),
      width=1600, height=1200, res=150)
  pheatmap(Xz, annotation_col = ann, show_rownames = TRUE, show_colnames = TRUE,
           main = paste0("Top ", min(topN, nrow(deg)), " DEGs — ", tag))
  dev.off()
}

dir.create("DEA/plots", showWarnings = FALSE)
heatmap_from_deg("DEA/DEG_ALAD_vs_Healthy_SIG_p0.05_lfc1.csv", "ALAD_vs_Healthy")
heatmap_from_deg("DEA/DEG_NLAD_vs_Healthy_SIG_p0.05_lfc1.csv", "NLAD_vs_Healthy")
heatmap_from_deg("DEA/DEG_ALAD_vs_NLAD_SIG_p0.05_lfc1.csv",    "ALAD_vs_NLAD")
