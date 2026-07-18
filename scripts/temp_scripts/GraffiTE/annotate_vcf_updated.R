#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(vcfR)

# Define the function to read and process the RepeatMasker output
read_rm_custom <- function(file) {
  rm_file <- readr::read_lines(file = file, skip = 3)
  rm_file <- lapply(rm_file, function(x) {
    str.res <- unlist(stringr::str_split(x, "\\s+"))
    str.res <- str.res[1:16]
    return(str.res)
  })
  rm_file <- tibble::as_tibble(do.call(rbind, rm_file))
  colnames(rm_file) <- c("sw_score", "perc_div", "perc_del",
                         "perc_insert", "qry_id", "qry_start", "qry_end", "qry_left", "strand",
                         "repeat_id", "matching_class", "in_repeat_start", "in_repeat_end", "in_repeat_left", "ID", "fragmts")
  qry_end <- qry_start <- NULL
  nrow_before_filtering <- nrow(rm_file)
  suppressWarnings(rm_file <- dplyr::mutate(rm_file,
                                            qry_start = as.integer(qry_start),
                                            qry_end = as.integer(qry_end), fragmts = as.integer(fragmts)
  ))
  rm_file <- dplyr::filter(rm_file, !is.na(qry_start), !is.na(qry_end))
  rm_file <- dplyr::mutate(rm_file, qry_width = as.integer(qry_end - qry_start + 1L))
  nrow_after_filtering <- nrow(rm_file)
  if ((nrow_before_filtering - nrow_after_filtering) > 0)
    message((nrow_before_filtering - nrow_after_filtering) +
              1, " out of ", nrow_before_filtering, " rows ~ ",
            round(((nrow_before_filtering - nrow_after_filtering) +
                     1) / nrow_before_filtering, 3), "% were removed from the imported RepeatMasker file, ",
            "because they contained 'NA' values in either 'qry_start' or 'qry_end'.")
  return(rm_file)
}

# Define paths
VCF <-  "/scratch/Users/olde5615/data/scale_up/bams_withID/sniffles2_variants.vcf"            # Path to your VCF file
REPMASK_ONECODE_OUT <- "/scratch/Users/olde5615/data/scale_up/bams_withID/repeatmasker_dirindels.fa.onecode.out"  # Path to RepeatMasker .out file
ANNOT_FILE <- "vcf_annotation_1"    # Path to the output annotation file

# Call the function to read the RepeatMasker file and assign to rm_tibble
rm_tibble <- read_rm_custom(REPMASK_ONECODE_OUT)

# Now rm_tibble contains your processed RepeatMasker data, and you can proceed with your analysis
print(rm_tibble)


# Read and process the RepeatMasker output
rep_mask <- read_rm_custom(REPMASK_ONECODE_OUT) %>%
  mutate(match_len = qry_end - qry_start) %>%
  arrange(qry_start) %>%
  group_by(qry_id) %>%
  summarise(
    start = first(qry_start), stop = max(qry_end),
    match_lengths = paste0(match_len, collapse = ","),
    fragmts = paste0(fragmts, collapse = ","),
    repeat_ids = paste0(repeat_id, collapse = ","),
    matching_classes = paste0(matching_class, collapse = ","),
    strands = paste0(strand, collapse = ","),
    RM_id = paste0(ID, collapse = ","),
    n_hits = n()
  )

# Load the VCF file and extract necessary information
vcf <- read.vcfR(VCF)
vcf_df <- tibble(
  CHROM = getCHROM(vcf),
  POS = getPOS(vcf),
  REF = getREF(vcf),
  ALT = getALT(vcf),
  qry_length = abs(str_length(getALT(vcf)) - str_length(getREF(vcf))),
  qry_id = getID(vcf)
)

# Merge VCF data with RepeatMasker annotation
annot <- left_join(vcf_df, rep_mask, by = "qry_id") %>%
  replace_na(list(matching_classes = "None",
                  repeat_ids = "None",
                  n_hits = 0,
                  fragmts = "0",
                  match_lengths = "0",
                  strands = "None",
                  RM_id = "None")) %>%
  select(-c(qry_length)) %>%
  arrange(CHROM, POS, qry_id) %>%
  select(CHROM, POS, qry_id, REF, ALT, n_hits, fragmts, match_lengths, repeat_ids, matching_classes, strands, RM_id)

# Write the annotation output file
write_tsv(annot, file = ANNOT_FILE, col_names = F)

# Print the column names for debugging purposes
print(colnames(annot))
