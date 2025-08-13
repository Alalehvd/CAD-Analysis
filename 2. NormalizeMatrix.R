# =======================
# Quantile-normalize a matrix CSV
# Input:  raw_probeset_matrix.csv  (rows=probesets, cols=samples)
# Output: normalized_dataset.csv (+ optional z-scored version)
# =======================

library(limma)

setwd("/Users/alalevd/Documents/GitHub/CAD-Analysis")

in_path  <- "raw_probeset_matrix.csv"
out_main <- "normalized_dataset.csv"
out_z    <- "normalized_dataset_zscore.csv"

# 1) Read matrix
X <- read.csv(in_path, row.names = 1, check.names = FALSE)

# 2) Coerce to numeric matrix
X <- as.matrix(data.frame(lapply(as.data.frame(X), as.numeric),
                          row.names = rownames(X), check.names = FALSE))

# 3) Decide whether to log2 (ONLY if clearly on linear scale)
#    expresso(medianpolish) already returns log2-scale by default.
rng <- range(X, na.rm = TRUE)
cat("Input range:", paste(signif(rng, 4), collapse = " .. "), "\n")

# Heuristic: if data look linear (very large values), do one log2; otherwise keep as-is
looks_linear <- (rng[2] > 30)  # MAS5 linear often >> 30; RMA/expresso(log2) typically <= ~20
if (looks_linear) {
  cat("Detected linear-scale input → applying single log2.\n")
  X_log2 <- log2(X + 1)
} else {
  cat("Detected log-scale input → skipping log2 (to avoid double-logging).\n")
  X_log2 <- X
}

# 4) Quantile normalization (columns = samples)
X_qn <- normalizeBetweenArrays(X_log2, method = "quantile")

# 5) Save normalized matrix
write.csv(X_qn, out_main)

# 6) (Optional) Z-score per gene (row-wise)
zscore_by_row <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[sdv == 0] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}
X_qn_z <- zscore_by_row(X_qn)
write.csv(X_qn_z, out_z)

cat("✅ Saved:\n -", out_main, "\n -", out_z, "\n")
