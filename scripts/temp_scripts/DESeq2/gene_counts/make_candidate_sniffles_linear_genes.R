library(stringr)
library(purrr)
library(tidyr)
library(dplyr)
library(DESeq2)
library(ggplot2)

inDir <- "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts_gene_id/perSample/"

# Get list of ALL featureCounts files
files_list <- list.files(path = inDir, pattern = "\\.txt$", full.names = FALSE)
head(files_list)

# Create list to store all sample data
df_list <- list()

# Loop through ALL files (no single-sample preprocessing needed)
for (i in seq_along(files_list)) {
  fn <- paste0(inDir, files_list[i])
  fc <- read.table(fn, header = TRUE)
  
  # Take first and last column (Geneid and counts)
  df_subset <- fc[, c(1, ncol(fc))]
  
  # Extract sample name from the last column name
  # This works for both the old (hg38.) and new (path) formats
  col_name <- as.character(colnames(df_subset)[2])
  
  # Get the basename (last part after the last slash)
  sample_name <- str_split(files_list[i], fixed(".featureCounts.txt"))[[1]][1]
  sample_name <- str_split(sample_name, fixed("_"))[[1]][1]  # Remove _2 
  
  colnames(df_subset) <- c("Geneid", sample_name)
  
  df_list[[i]] <- df_subset
}
head(df_list,n=5)
head(sample_name)
# Combine all samples into one count matrix
counts_df <- purrr::reduce(df_list, full_join, by = "Geneid")
head(counts_df,n=5)
write.csv(counts_df, "RNAseq_counts_allpeople_linear_genes.csv", row.names = FALSE)
