# ===== PCA with metadata, colors, clean labels =====
# Inputs: normalized_dataset.csv (matrix), samples_metadata.csv with columns:
#   GSM, file, group, [optional] batch
# Outputs: PCA_by_group_repel.png/pdf, PCA_by_group_SHORTLABELS.png/pdf, sample_label_map.csv

setwd("/Users/alalevd/Documents/GitHub/CAD-Analysis")

# ---- Setup ----
in_expr <- "normalized_dataset.csv"   
in_meta <- "samples_metadata.csv"
out_dir <- "QC_normalized2"
dir.create(out_dir, showWarnings = FALSE)

# Packages
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly=TRUE)) install.packages("ggrepel")
library(ggplot2)
library(ggrepel)

# ---- Load expression ----
X <- read.csv(in_expr, row.names = 1, check.names = FALSE)
X <- as.matrix(data.frame(lapply(as.data.frame(X), as.numeric), row.names=rownames(X)))
X <- X[rowSums(is.na(X)) < ncol(X), , drop=FALSE]

# ---- Load metadata ----
meta <- read.csv(in_meta, stringsAsFactors = FALSE)
needed <- c("GSM","file","group")
if (!all(needed %in% names(meta))) {
  stop("samples_metadata.csv must have columns: GSM, file, group")
}

# Align metadata order to expression columns
cn <- colnames(X)
idx_gsm <- match(cn, meta$GSM)

# Fallback match by file name if needed
need_file_match <- is.na(idx_gsm)
if (any(need_file_match)) {
  idx_file <- match(cn[need_file_match], meta$file)
  idx_gsm[need_file_match] <- idx_file
}

# Fallback match by extracting GSM from colnames
need_regex <- is.na(idx_gsm)
if (any(need_regex)) {
  cn_gsm <- sub(".*\\b(GSM\\d+)\\b.*", "\\1", cn[need_regex])
  idx_rgx <- match(cn_gsm, meta$GSM)
  idx_gsm[need_regex] <- idx_rgx
}

# Final check
if (any(is.na(idx_gsm))) {
  bad <- cn[is.na(idx_gsm)]
  stop(paste0("Could not map these columns to metadata: ", paste(bad, collapse=", ")))
}

meta_aligned <- meta[idx_gsm, , drop=FALSE]
colnames(X) <- meta_aligned$GSM

# ---- Short labels (S01, S02, …) ----
short_labels <- sprintf("S%02d", seq_len(ncol(X)))
label_map <- data.frame(
  short_label = short_labels,
  GSM = meta_aligned$GSM,
  file = meta_aligned$file,
  group = meta_aligned$group,
  stringsAsFactors = FALSE
)
if ("batch" %in% names(meta_aligned)) {
  label_map$batch <- meta_aligned$batch
}
write.csv(label_map, file.path(out_dir, "sample_label_map.csv"), row.names = FALSE)

# ---- PCA ----
pc <- prcomp(t(X), scale.=TRUE)
pvar <- summary(pc)$importance[2, 1:2] * 100
pca_df <- data.frame(
  PC1 = pc$x[,1],
  PC2 = pc$x[,2],
  GSM = colnames(X),
  group = factor(meta_aligned$group),
  stringsAsFactors = FALSE
)
if ("batch" %in% names(meta_aligned)) {
  pca_df$batch <- factor(meta_aligned$batch)
}
pca_df$short <- short_labels

# ---- Plot 1: GSM labels ----
p1 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.4, aes(shape = if ("batch" %in% names(pca_df)) batch else NULL)) +
  geom_text_repel(aes(label = GSM), size = 3, max.overlaps = Inf,
                  box.padding = 0.5, point.padding = 0.2) +
  stat_ellipse(aes(color = group), level = 0.68, linetype = 2) +
  labs(title = "PCA (normalized) by group",
       x = paste0("PC1 (", round(pvar[1], 1), "%)"),
       y = paste0("PC2 (", round(pvar[2], 1), "%)")) +
  theme_minimal()

ggsave(file.path(out_dir, "PCA_by_group_repel.pdf"), p1, width = 10, height = 8)
ggsave(file.path(out_dir, "PCA_by_group_repel.png"), p1, width = 10, height = 8, dpi = 300)

# ---- Plot 2: Short labels ----
p2 <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.8, aes(shape = if ("batch" %in% names(pca_df)) batch else NULL)) +
  geom_text_repel(aes(label = short), size = 3, max.overlaps = Inf,
                  box.padding = 0.5, point.padding = 0.2) +
  stat_ellipse(aes(color = group), level = 0.68, linetype = 2) +
  labs(title = "PCA (normalized) by group — short labels",
       x = paste0("PC1 (", round(pvar[1], 1), "%)"),
       y = paste0("PC2 (", round(pvar[2], 1), "%)")) +
  theme_minimal()

ggsave(file.path(out_dir, "PCA_by_group_SHORTLABELS.pdf"), p2, width = 10, height = 8)
ggsave(file.path(out_dir, "PCA_by_group_SHORTLABELS.png"), p2, width = 10, height = 8, dpi = 300)

cat("✅ Saved PCA plots and sample_label_map.csv in", out_dir, "\n")
