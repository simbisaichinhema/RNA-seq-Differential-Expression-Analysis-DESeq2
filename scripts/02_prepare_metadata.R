############################################################
# 02_prepare_metadata.R
# Purpose: Parse GEO series matrix metadata to create a
#          structured sample metadata table for downstream
#          DESeq2 analysis. Filters to the comparison of
#          interest: "Healthy control" vs "No vaccination,
#          No prior infection" (unvaccinated infected).
#
# Required packages: None (base R only)
#
# Inputs:
#   data/raw/GSE201530_series_matrix.txt
#
# Outputs:
#   data/metadata/Metadata.csv
#
# Expected runtime: < 1 minute
#
# Example execution:
#   Rscript scripts/02_prepare_metadata.R
#
# Error handling: Stops with informative message if metadata
#   parsing fails or samples are unmatched.
############################################################

cat("========================================\n")
cat(" Step 2: Prepare Metadata\n")
cat("========================================\n\n")

############################################################
# Read GEO Series Matrix
############################################################

cat("Reading GEO series matrix...\n")

series_matrix_file <- "data/raw/GSE201530_series_matrix.txt"

if (!file.exists(series_matrix_file)) {
  stop(paste0("Series matrix file not found: ", series_matrix_file))
}

geo <- readLines(series_matrix_file)

cat("Done.\n\n")

############################################################
# Extract GSM IDs
############################################################

gsm <- sub("!Sample_geo_accession = ", "", geo[grep("^!Sample_geo_accession", geo)])
gsm <- strsplit(gsm, "\t")[[1]]
gsm <- gsub('"', "", gsm)

cat("Found", length(gsm), "samples.\n")

############################################################
# Extract Sample Titles
############################################################

titles <- sub("!Sample_title = ", "", geo[grep("^!Sample_title", geo)])
sample_names <- strsplit(titles, "\t")[[1]]
sample_names <- gsub('"', "", sample_names)

############################################################
# Extract Characteristics
############################################################

char_lines <- geo[grep("^!Sample_characteristics_ch1", geo)]

extract_metadata <- function(line) {
  values <- strsplit(line, "\t")[[1]][-1]
  values <- gsub('"', "", values)
  field <- sub(":.*", "", values[1])
  field <- trimws(field)
  values <- sub("^[^:]+: ", "", values)
  list(field = field, values = values)
}

metadata_list <- lapply(char_lines, extract_metadata)
names(metadata_list) <- sapply(metadata_list, function(x) x$field)

############################################################
# Build Metadata Table
############################################################

n_samples <- length(gsm)
for (nm in names(metadata_list)) {
  n <- length(metadata_list[[nm]]$values)
  if (n < n_samples) {
    metadata_list[[nm]]$values <- c(
      metadata_list[[nm]]$values,
      rep(NA_character_, n_samples - n)
    )
  }
}
rm(nm, n, n_samples)

metadata <- data.frame(
  GSM = gsm,
  Sample = sample_names,
  Gender = metadata_list[["gender"]]$values,
  Age = as.numeric(metadata_list[["age"]]$values),
  Group = metadata_list[["group (by covid-19 vaccination, prior infection)"]]$values,
  Omicron_Lineage = metadata_list[["omicron sublineage"]]$values,
  Days_After_PCR = metadata_list[["days after positive pcr results"]]$values,
  Disease_State = metadata_list[["disease state"]]$values,
  Country = metadata_list[["geographical location"]]$values,
  Cell_Type = metadata_list[["cell type"]]$values,
  stringsAsFactors = FALSE
)

############################################################
# Filter to Comparison of Interest
############################################################

cat("\nFiltering to comparison groups:\n")

selected_groups <- c("Healthy control", "No vaccination, No prior infection")

metadata_filtered <- metadata[metadata$Group %in% selected_groups, ]

metadata_filtered$Group <- factor(
  metadata_filtered$Group,
  levels = c("Healthy control", "No vaccination, No prior infection")
)

cat("Healthy control:", nrow(metadata_filtered[metadata_filtered$Group == "Healthy control", ]), "samples\n")
cat("Unvaccinated infected:", nrow(metadata_filtered[metadata_filtered$Group == "No vaccination, No prior infection", ]), "samples\n")

############################################################
# Save Metadata
############################################################

output_file <- "data/metadata/Metadata.csv"
write.csv(metadata_filtered, output_file, row.names = FALSE)

cat("\nMetadata saved to:", output_file, "\n")
cat("Total samples:", nrow(metadata_filtered), "\n\n")

############################################################
# Verification
############################################################

if (nrow(metadata_filtered) == 0) {
  stop("No samples matched the selected groups. Check metadata column names.")
}

if (any(is.na(metadata_filtered$Group))) {
  warning("Some samples have missing group assignments.")
}

cat("========================================\n")
cat(" Step 2 Complete\n")
cat("========================================\n")
