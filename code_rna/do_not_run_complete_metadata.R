#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript complete_sample_metadata.R <metadata_csv> <fastq_dir>")
}

metadata_file <- args[1]
fastq_dir <- args[2]

library(dplyr)
library(readr)
library(stringr)

# Read metadata
meta <- read_csv(metadata_file, show_col_types = FALSE)

# Check required columns
required_cols <- c("sample_id", "condition", "file_name_1")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("❌ Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# Check if file_name_2 exists (optional for paired-end)
has_file_2 <- "file_name_2" %in% names(meta)

# Function to search for file in directory
find_fastq_path <- function(filename) {
  if (is.na(filename) || filename == "") return(NA_character_)
  found <- list.files(fastq_dir, pattern = paste0("^", filename, "$"), full.names = TRUE, recursive = TRUE)
  if (length(found) == 0) stop(paste("❌ File not found in FASTQ directory:", filename))
  return(found[1])
}

# Resolve paths
meta <- meta %>%
  mutate(
    file_path_1 = sapply(file_name_1, find_fastq_path),
    file_path_2 = if (has_file_2) sapply(file_name_2, find_fastq_path) else NA_character_
  )

# Add replicate per condition
meta <- meta %>%
  group_by(condition) %>%
  mutate(replicate = paste0("rep", row_number())) %>%
  ungroup()

# Save updated metadata
write_csv(meta, metadata_file)
