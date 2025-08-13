library(affy)

# Read the data
cel_dir <- "/Users/alalevd/Documents/GitHub/CAD-Analysis/Raw"
abatch  <- ReadAffy(celfile.path = cel_dir)

# Summarization 
eset_rawsum <- expresso(
  abatch,
  bgcorrect.method = "none",
  normalize        = FALSE,
  pmcorrect.method = "pmonly",
  summary.method   = "medianpolish"  # یا "average" اگر medianpolish نخواستی
)

# Create the matrix
mat_rawsum <- exprs(eset_rawsum)  # rows = probeset IDs

# Save as CSV
proj <- "/Users/alalevd/Documents/GitHub/CAD-Analysis"
dir.create(proj, recursive = TRUE, showWarnings = FALSE)

write.csv(mat_rawsum, file.path(proj, "raw_probeset_matrix.csv"))
write.csv(log2(mat_rawsum + 1), file.path(proj, "log2_raw_probeset_matrix.csv"))
