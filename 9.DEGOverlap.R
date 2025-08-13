library(UpSetR)

dir.create("DEA2/plots", showWarnings = FALSE)

# === 3) Overlap (UpSet) of significant gene sets ===
deg_ALAD <- read.csv("DEA2/DEG_ALAD_vs_Healthy_SIG_p0.05_lfc1.csv")
deg_NLAD <- read.csv("DEA2/DEG_NLAD_vs_Healthy_SIG_p0.05_lfc1.csv")
deg_ALAD_NLAD <- read.csv("DEA2/DEG_ALAD_vs_NLAD_SIG_p0.05_lfc1.csv")

# Extract gene names
genes_ALAD <- rownames(deg_ALAD)
genes_NLAD <- rownames(deg_NLAD)
genes_ALAD_NLAD <- rownames(deg_ALAD_NLAD)

# Put them into a list
gene_lists <- list(
  ALAD_vs_Healthy = genes_ALAD,
  NLAD_vs_Healthy = genes_NLAD,
  ALAD_vs_NLAD    = genes_ALAD_NLAD
)

# Save overlap as CSV in DEA/plots
all_genes <- unique(unlist(gene_lists))
overlap_df <- data.frame(
  Gene = all_genes,
  ALAD_vs_Healthy = all_genes %in% genes_ALAD,
  NLAD_vs_Healthy = all_genes %in% genes_NLAD,
  ALAD_vs_NLAD    = all_genes %in% genes_ALAD_NLAD
)
write.csv(overlap_df, file.path("DEA2/plots", "overlap_matrix.csv"), row.names = FALSE)

# Draw UpSet plot in DEA/plots
png(file.path("DEA2/plots", "upset_overlap.png"), width = 1200, height = 900, res = 150)
upset(fromList(gene_lists),
      order.by = "freq",
      mainbar.y.label = "Gene Overlap",
      sets.x.label = "Genes per Contrast")
dev.off()
