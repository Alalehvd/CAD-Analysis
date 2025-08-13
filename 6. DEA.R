library(limma)

# ====== PARAMETERS ======
PADJ_METHOD  <- "BH"   # "BH", "bonferroni", "holm", "BY", ...
PADJ_CUTOFF  <- 0.05   # adjusted p-value threshold
LOGFC_CUTOFF <- 1      # absolute log2 fold-change threshold
TOP_N        <- Inf    # how many rows to return from topTable (Inf = all)
OUT_DIR      <- "DEA2"  # folder for all output
dir.create(OUT_DIR, showWarnings = FALSE)
# ========================

# ---- Load ----
expr <- read.csv("normalized_dataset.csv", row.names = 1, check.names = FALSE)
expr <- as.matrix(data.frame(lapply(as.data.frame(expr), as.numeric),
                             row.names = rownames(expr), check.names = FALSE))

meta <- read.csv("samples_metadata_grouped.csv", stringsAsFactors = FALSE)

# ---- 1) Rename columns to GSM first (from filenames like 'GSMxxxxx_....CEL.gz') ----
gsm_ids <- sub("^(GSM\\d+).*", "\\1", colnames(expr))
colnames(expr) <- gsm_ids

# ---- 2) Match metadata to expression by GSM ----
idx <- match(colnames(expr), meta$GSM)
if (any(is.na(idx))) {
  stop("Some GSM IDs in expression are not found in metadata: ",
       paste(colnames(expr)[is.na(idx)], collapse = ", "))
}
meta <- meta[idx, , drop = FALSE]

# ---- 3) Factor group & sanity checks ----
meta$group <- factor(meta$group, levels = c("Healthy","NLAD","ALAD"))
stopifnot(ncol(expr) == nrow(meta))

# ---- 4) Design & contrasts ----
design <- model.matrix(~ 0 + meta$group)
colnames(design) <- levels(meta$group)
rownames(design) <- colnames(expr)

contrast.matrix <- makeContrasts(
  ALAD_vs_Healthy = ALAD - Healthy,
  NLAD_vs_Healthy = NLAD - Healthy,
  ALAD_vs_NLAD    = ALAD - NLAD,
  levels = design
)

# ---- 5) Fit model ----
fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# ---- 6) Save results helper ----
save_results <- function(fit2, coef_name, base_name) {
  tt  <- topTable(fit2, coef = coef_name, number = TOP_N, adjust.method = PADJ_METHOD)
  write.csv(tt, file.path(OUT_DIR, paste0("DEG_", base_name, "_ALL.csv")))
  
  sig <- subset(tt, adj.P.Val < PADJ_CUTOFF & abs(logFC) > LOGFC_CUTOFF)
  write.csv(sig, file.path(OUT_DIR, paste0("DEG_", base_name,
                                           "_SIG_p", PADJ_CUTOFF,
                                           "_lfc", LOGFC_CUTOFF, ".csv")))
  cat(base_name, "→", nrow(sig), "significant genes at adj.p <", PADJ_CUTOFF,
      "and |logFC| >", LOGFC_CUTOFF, "\n")
  invisible(list(all = tt, sig = sig))
}

# ---- 7) Run contrasts ----
res_ALAD_Healthy <- save_results(fit2, "ALAD_vs_Healthy", "ALAD_vs_Healthy")
res_NLAD_Healthy <- save_results(fit2, "NLAD_vs_Healthy", "NLAD_vs_Healthy")
res_ALAD_NLAD    <- save_results(fit2, "ALAD_vs_NLAD",    "ALAD_vs_NLAD")

cat("\n✅ All DEA results saved in:", normalizePath(OUT_DIR), "\n")
