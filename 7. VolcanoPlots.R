library(ggplot2)

dir.create("DEA/plots", showWarnings = FALSE)

make_volcano <- function(file, label, padj_cut=0.05, lfc_cut=1) {
  df <- read.csv(file, row.names = 1, check.names = FALSE)
  df$negLog10FDR <- -log10(pmax(df$adj.P.Val, .Machine$double.xmin))
  
  # Classify for color
  df$Status <- "Not Sig"
  df$Status[df$adj.P.Val < padj_cut & df$logFC >  lfc_cut] <- "Up"
  df$Status[df$adj.P.Val < padj_cut & df$logFC < -lfc_cut] <- "Down"
  df$Status <- factor(df$Status, levels = c("Up", "Down", "Not Sig"))
  
  p <- ggplot(df, aes(x=logFC, y=negLog10FDR, color=Status)) +
    geom_point(alpha=0.6, size=1.5) +
    scale_color_manual(values = c("Up"="red", "Down"="blue", "Not Sig"="grey")) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed", color="black") +
    geom_hline(yintercept=-log10(padj_cut), linetype="dashed", color="black") +
    labs(title=paste("Volcano Plot —", label),
         x="log2 Fold Change",
         y="-log10(Adjusted P-value)") +
    theme_minimal(base_size = 14) +
    theme(legend.position="top")
  
  ggsave(filename=sub("\\.csv$", "_volcano.png", 
                      file.path("DEA/plots", basename(file))),
         plot=p, width=8, height=6, dpi=300)
}

make_volcano("DEA/DEG_ALAD_vs_Healthy_ALL.csv",   "ALAD vs Healthy")
make_volcano("DEA/DEG_NLAD_vs_Healthy_ALL.csv",   "NLAD vs Healthy")
make_volcano("DEA/DEG_ALAD_vs_NLAD_ALL.csv",      "ALAD vs NLAD")
