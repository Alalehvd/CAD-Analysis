library(affy)
library(canine2cdf) 

cel_dir <- "/Users/alalehvd/Documents/GitHub/CAD-Analysis/Raw"

# Read all CEL(.gz)
data_affy <- ReadAffy(celfile.path = cel_dir)

# RMA normalization (log2)
eset <- rma(data_affy, verbose = TRUE)

# Save expression matrix
expr <- exprs(eset)
write.csv(expr, "/Users/alalehvd/Documents/GitHub/CAD-Analysis/expression_rma_probes.csv")
cat("✅ Saved: /Users/alalehvd/Documents/GitHub/CAD-Analysis/expression_rma_probes.csv\n")

# Quick sanity checks (optional)
dim(expr)               # rows=probes, cols=samples
colnames(expr)[1:5]     # first few sample IDs
rownames(expr)[1:5]     # first few probe IDs


# ---- Build sample metadata from column names and save ----

# 1) Build metadata from real CEL filenames
files <- list.files(cel_dir, pattern = "\\.CEL(\\.gz)?$", full.names = FALSE)

# GSM from filename
gsm_from_file <- sub("^(GSM\\d+).*", "\\1", files)

# Group: check Nonlesional BEFORE Lesional
group_from_file <- ifelse(grepl("Nonlesional", files, ignore.case = TRUE), "NLAD",
                          ifelse(grepl("Lesional", files, ignore.case = TRUE), "ALAD",
                                 ifelse(grepl("Normal|Healthy", files, ignore.case = TRUE), "Healthy", "Unknown")))

meta_all <- data.frame(
  GSM  = gsm_from_file,
  file = files,
  group = group_from_file,
  stringsAsFactors = FALSE
)

# 2) Load the GSM-named expression matrix and align
expr <- read.csv("/Users/alalehvd/Documents/GitHub/CAD-Analysis/expression_rma_probes_GSMcols.csv",
                 row.names = 1, check.names = FALSE)

# Keep only metadata rows that exist in expr, and in the same order
stopifnot(all(colnames(expr) %in% meta_all$GSM))
meta <- meta_all[match(colnames(expr), meta_all$GSM), ]

# 3) Sanity check
print(table(meta$group))  # should be 13 Healthy, 13 NLAD, 13 ALAD

# 4) Save corrected metadata
write.csv(meta, "/Users/alalehvd/Documents/GitHub/CAD-Analysis/samples_metadata.csv", row.names = FALSE)

# 5) (Optional) Save expression again to be safe
write.csv(expr, "/Users/alalehvd/Documents/GitHub/CAD-Analysis/expression_rma_probes_GSMcols.csv")
cat("✅ Metadata fixed and saved.\n")

