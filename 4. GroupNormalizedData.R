# ---- Setup ----
expr_file <- "normalized_dataset.csv"
meta_file <- "samples_metadata.csv"

# Read expression data
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr <- as.matrix(data.frame(lapply(as.data.frame(expr), as.numeric), 
                             row.names = rownames(expr)))

# Read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Ensure metadata order matches expression columns
# First try GSM match
idx <- match(colnames(expr), meta$GSM)

# If NA, try matching by filename
if (any(is.na(idx))) {
  idx_file <- match(colnames(expr), meta$file)
  idx[is.na(idx)] <- idx_file[is.na(idx)]
}

# Reorder metadata to match expression
meta <- meta[idx, ]
stopifnot(all(!is.na(meta$group)))  # make sure every sample got a group

# Convert group to factor (set desired order)
meta$group <- factor(meta$group, levels = c("Healthy", "NLAD", "ALAD"))

# Quick check
table(meta$group)

# Save grouped expression + metadata
write.csv(meta, "samples_metadata_grouped.csv", row.names = FALSE)
write.csv(expr, "normalized_dataset_grouped.csv")

cat("✅ Grouping complete. Metadata and expression saved.\n")
