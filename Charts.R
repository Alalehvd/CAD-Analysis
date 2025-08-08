#########################################
# QC & Exploratory Analysis for CAD data
#########################################

library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# 1. Boxplot of expression values per sample
pdf(file.path(out_dir, "QC_boxplot_expression.pdf"), width = 10, height = 5)
boxplot(expr,
        main = "Boxplot of Expression Values",
        xlab = "Samples", ylab = "Expression (log2)",
        las = 2, outline = FALSE,
        col = brewer.pal(8, "Set3")[as.numeric(meta$group)])
legend("topright", legend = levels(meta$group),
       fill = brewer.pal(8, "Set3")[1:length(levels(meta$group))])
dev.off()

# 2. Density plot of expression distributions
pdf(file.path(out_dir, "QC_density_expression.pdf"), width = 7, height = 5)
plot(density(expr[,1]), col = brewer.pal(8, "Set3")[as.numeric(meta$group[1])],
     main = "Density Plot of Expression Values", xlab = "Expression (log2)")
for (i in 2:ncol(expr)) {
  lines(density(expr[,i]),
        col = brewer.pal(8, "Set3")[as.numeric(meta$group[i])])
}
legend("topright", legend = levels(meta$group),
       fill = brewer.pal(8, "Set3")[1:length(levels(meta$group))])
dev.off()

# 3. PCA plot (top 500 most variable genes)
top_var_genes <- order(apply(expr, 1, var), decreasing = TRUE)[1:500]
pca <- prcomp(t(expr[top_var_genes, ]), scale. = TRUE)

pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     Group = meta$group)

percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))[1:2]

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_minimal() +
  ggtitle("PCA - Top 500 Variable Genes") +
  scale_color_brewer(palette = "Set1")
ggsave(file.path(out_dir, "QC_PCA_plot.png"), plot = p, width = 6, height = 5, dpi = 150)

# 4. Sample distance heatmap
sample_dist <- dist(t(expr[top_var_genes, ]))
sample_dist_mat <- as.matrix(sample_dist)
rownames(sample_dist_mat) <- colnames(expr)
colnames(sample_dist_mat) <- colnames(expr)

pdf(file.path(out_dir, "QC_sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(sample_dist_mat,
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Distance Heatmap")
dev.off()

#########################################
# 1. Mean–variance plot (limma)
#########################################
png(file.path(out_dir, "mean_variance_plot.png"), width = 1200, height = 900, res = 150)
plotSA(fit2, main = "Mean–Variance Trend (limma)")
dev.off()

#########################################
# 2. TopTable plot - ALAD vs Healthy
#########################################
top_ALAD <- topTable(fit2, coef = "ALADvsHealthy", number = 20, sort.by = "P")
png(file.path(out_dir, "top20_ALADvsHealthy_bar.png"), width = 1200, height = 900, res = 150)
barplot(-log10(top_ALAD$adj.P.Val),
        names.arg = top_ALAD$GeneSymbol,
        las = 2, cex.names = 0.7,
        col = "steelblue",
        main = "Top 20 DEGs (ALAD vs Healthy)",
        ylab = "-log10 Adjusted P-value")
dev.off()

#########################################
# 3. UpSet plot for overlap
#########################################
library(UpSetR)

# Significant gene sets
sig_ALAD_genes <- subset(deg_ALAD, adj.P.Val < 0.05 & abs(logFC) > 1)$GeneSymbol
sig_NLAD_genes <- subset(deg_NLAD, adj.P.Val < 0.05 & abs(logFC) > 1)$GeneSymbol

# Combine into list
gene_lists <- list(ALAD_vs_Healthy = sig_ALAD_genes,
                   NLAD_vs_Healthy = sig_NLAD_genes)

# Remove empty sets
gene_lists <- Filter(function(x) length(x) > 0, gene_lists)

# Only plot if at least 2 sets have data
if (length(gene_lists) >= 2) {
  png(file.path(out_dir, "upset_overlap.png"), width = 1200, height = 900, res = 150)
  upset(fromList(gene_lists),
        order.by = "freq", 
        mainbar.y.label = "Gene Overlap",
        sets.x.label = "Genes per Contrast")
  dev.off()
  cat("✅ UpSet plot saved to:", file.path(out_dir, "upset_overlap.png"), "\n")
} else {
  cat("⚠️ Not enough non-empty sets to create UpSet plot. Skipped.\n")
}

#########################################
# 4. Venn diagram
#########################################
library(limma)  # vennDiagram function

# Create logical matrix for Venn
venn_input <- cbind(
  ALAD = rownames(expr) %in% sig_ALAD_genes,
  NLAD = rownames(expr) %in% sig_NLAD_genes
)

png(file.path(out_dir, "venn_overlap.png"), width = 900, height = 900, res = 150)
vennDiagram(venn_input, names = c("ALAD", "NLAD"),
            circle.col = c("skyblue", "pink"))
dev.off()
