library(stringr)
library(purrr)
library(tidyr)
library(dplyr)
library(DESeq2)
library(ggplot2)

inDir <- "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/perSample/"

# Get list of ALL featureCounts files
files_list <- list.files(path = inDir, pattern = "\\.txt$", full.names = FALSE)

# Create list to store all sample data
df_list <- list()
# Loop through ALL files (no single-sample preprocessing needed)
for (i in seq_along(files_list)) {   
  fn <- paste0(inDir, files_list[i])
  fc <- read.table(fn, header = TRUE)
  
  # Take first and last column (Geneid and counts)
  df_subset <- fc[, c(1, ncol(fc))]
  
  # Extract sample name from column name (removing hg38. and .merged.bam)
  newname <- str_split(as.character(colnames(df_subset)[2]), pattern = "hg38.")[[1]][2]
  newname <- str_split(newname, pattern = fixed("."))[[1]][1]
  colnames(df_subset) <- c("Geneid", newname)
  
  df_list[[i]] <- df_subset 
}

# Combine all samples into one count matrix
counts_df <- purrr::reduce(df_list, full_join, by = "Geneid")

write.csv(counts_df, "RNAseq_counts_allpeople.csv")

