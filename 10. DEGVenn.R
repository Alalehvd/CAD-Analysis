# ===== Venn diagram for significant DEGs (3 contrasts) =====
# Inputs (already created by your DEA step):
#   DEA/DEG_ALAD_vs_Healthy_SIG_p0.05_lfc1.csv
#   DEA/DEG_NLAD_vs_Healthy_SIG_p0.05_lfc1.csv
#   DEA/DEG_ALAD_vs_NLAD_SIG_p0.05_lfc1.csv
# Output:
#   DEA/plots/venn_sig_DEGs.png

if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)

dir.create("DEA2/plots", showWarnings = FALSE)

# Load significant gene lists (row names are genes)
read_sig <- function(path) {
  if (!file.exists(path)) return(character())
  x <- read.csv(path, row.names = 1, check.names = FALSE)
  rownames(x)
}

genes_AH <- read_sig("DEA2/DEG_ALAD_vs_Healthy_SIG_p0.05_lfc1.csv")
genes_NH <- read_sig("DEA2/DEG_NLAD_vs_Healthy_SIG_p0.05_lfc1.csv")
genes_AN <- read_sig("DEA2/DEG_ALAD_vs_NLAD_SIG_p0.05_lfc1.csv")

sets <- list("ALAD vs Healthy" = genes_AH,
             "NLAD vs Healthy" = genes_NH,
             "ALAD vs NLAD"    = genes_AN)

# If all sets are empty, skip gracefully
if (all(lengths(sets) == 0)) {
  message("No significant DEGs to plot in Venn.")
} else {
  # VennDiagram requires at least 2 non-empty sets
  nnz <- sum(lengths(sets) > 0)
  if (nnz < 2) {
    message("Need at least 2 non-empty sets for a Venn diagram.")
  } else {
    # Draw and save
    png("DEA2/plots/venn_sig_DEGs.png", width = 1600, height = 1200, res = 150)
    grid.newpage()
    v <- venn.diagram(
      x = sets,
      filename = NULL,                 # draw to current device
      fill = c("#F94144", "#277DA1", "#90BE6D"),  # colors
      alpha = 0.5,
      lwd = 2,
      cex = 1.4,                       # numbers size
      cat.cex = 1.4,                   # labels size
      cat.pos = 0,
      cat.dist = 0.05,
      margin = 0.08
    )
    grid.draw(v)
    dev.off()
    message("✅ Saved: DEA2/plots/venn_sig_DEGs.png")
  }
}

