# ===== QC for normalized microarray matrix =====
# Input: normalized_dataset.csv (rows=genes/probes, cols=samples, log2-scale)
# Outputs: PDFs/PNGs with QC plots

setwd("/Users/alalevd/Documents/GitHub/CAD-Analysis")

# ---- Setup ----
in_file <- "normalized_dataset.csv"           # change path if needed
out_dir <- "QC_normalized2"
dir.create(out_dir, showWarnings = FALSE)

# Packages
pkg_ok <- function(p) { if (!requireNamespace(p, quietly=TRUE)) install.packages(p); library(p, character.only=TRUE) }
pkg_ok("pheatmap")

# ---- Load ----
X <- read.csv(in_file, row.names = 1, check.names = FALSE)
# coerce to numeric (protect against accidental character columns)
X <- as.matrix(data.frame(lapply(as.data.frame(X), as.numeric), row.names=rownames(X)))
# drop rows with all NA
X <- X[rowSums(is.na(X)) < ncol(X), , drop=FALSE]

# ---- Rename columns from file names to GSM IDs ----
colnames(X) <- sub("^(GSM\\d+).*", "\\1", colnames(X))

# ---- Quick checks ----
cat("Dims:", paste(dim(X), collapse=" x "), "\n")
cat("Range:", range(X, na.rm=TRUE), "\n")   

# ---- Boxplot (per sample) ----
pdf(file.path(out_dir, "boxplot_samples.pdf"), width=10, height=6)
par(mar=c(8,4,2,1))
boxplot(X, las=2, outline=FALSE, main="Boxplot per sample (normalized)")
dev.off()

# ---- Density overlay ----
pdf(file.path(out_dir, "density_samples.pdf"), width=10, height=6)
plot(density(na.omit(X[,1])), main="Density per sample (normalized)", xlab="Expression", lwd=2)
for (j in 2:ncol(X)) lines(density(na.omit(X[,j])), lwd=1)
legend("topright", legend=c("All samples"), bty="n")
dev.off()

# ---- PCA ----
X_t <- t(X)  # samples as rows
pc <- prcomp(X_t, scale.=TRUE)
pca_df <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], sample=colnames(X))

pdf(file.path(out_dir, "PCA_samples.pdf"), width=12, height=9)  # bigger canvas
plot(pca_df$PC1, pca_df$PC2, pch=19,
     xlab=paste0("PC1 (", round(100*summary(pc)$importance[2,1],1), "%)"),
     ylab=paste0("PC2 (", round(100*summary(pc)$importance[2,2],1), "%)"),
     main="PCA (normalized)")
text(pca_df$PC1, pca_df$PC2, labels=pca_df$sample, pos=3, cex=0.6)  # smaller
dev.off()


# ---- Sample correlation heatmap ----
cors <- cor(X, use="pairwise.complete.obs")
pheatmap::pheatmap(cors, main="Sample correlation (normalized)",
                   filename=file.path(out_dir, "cor_heatmap.png"), width=8, height=8)

# ---- Hierarchical clustering dendrogram ----
pdf(file.path(out_dir, "hclust_samples.pdf"), width=7, height=6)
plot(hclust(as.dist(1 - cors), method="average"), main="Sample clustering (1 - correlation)")
dev.off()

cat("✅ QC plots saved in:", normalizePath(out_dir), "\n")
