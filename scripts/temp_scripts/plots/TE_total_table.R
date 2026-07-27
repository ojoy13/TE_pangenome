# Load required libraries
library(tidyverse)
library(vcfR)

# Create output directory
output_dir <- "INS_DEL_plot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created directory:", output_dir, "\n")
}

# Read the VCF file
vcf_file <- "/scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/pangenome_21_dedup.vcf"
vcf <- read.vcfR(vcf_file)

# Extract INFO fields
info_df <- as.data.frame(vcf@fix[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")])

# Extract sample names (column headers after INFO)
sample_names <- colnames(vcf@gt)[-1]  # Remove the FORMAT column
n_samples <- length(sample_names)
cat("Number of samples:", n_samples, "\n")

# Simple function to extract values from INFO using regex
extract_info_value <- function(info_string, key) {
  pattern <- paste0(key, "=([^;]+)")
  match <- regexpr(pattern, info_string)
  if (match[1] > 0) {
    start <- match[1] + nchar(key) + 1
    end <- start + attr(match, "match.length") - nchar(key) - 2
    return(substr(info_string, start, end))
  } else {
    return(NA)
  }
}

# Extract SVTYPE and repeat_ids
cat("Parsing VCF INFO fields...\n")
variants_df <- data.frame(
  CHROM = info_df$CHROM,
  POS = info_df$POS,
  ID = info_df$ID,
  SVTYPE = sapply(info_df$INFO, function(x) extract_info_value(x, "SVTYPE")),
  repeat_ids = sapply(info_df$INFO, function(x) extract_info_value(x, "repeat_ids")),
  stringsAsFactors = FALSE
)

# Extract genotype information from VCF
cat("Extracting genotype information...\n")
gt_matrix <- vcf@gt[, -1]  # Remove FORMAT column
colnames(gt_matrix) <- sample_names

# Parse genotypes - count ALT alleles (0, 1, or 2)
parse_genotype <- function(gt_string) {
  if (is.na(gt_string) || gt_string == "." || gt_string == "./." || gt_string == ".|.") {
    return(NA)  # Missing genotype
  } else if (gt_string == "0/0" || gt_string == "0|0") {
    return(0)  # Homozygous reference
  } else if (gt_string == "1/1" || gt_string == "1|1" || 
             gt_string == "2/2" || gt_string == "2|2" ||
             gt_string == "3/3" || gt_string == "3|3" ||
             gt_string == "4/4" || gt_string == "4|4") {
    return(2)  # Homozygous alternate
  } else if (grepl("0/1|1/0|0\\|1|1\\|0", gt_string)) {
    return(1)  # Heterozygous
  } else {
    # For complex genotypes, count non-zero alleles
    alleles <- unlist(strsplit(gt_string, "[/|]"))
    alt_count <- sum(alleles != "0")
    return(alt_count)
  }
}

# Apply genotype parsing to all samples
for (sample in sample_names) {
  variants_df[[paste0(sample, "_ALT")]] <- sapply(gt_matrix[, sample], parse_genotype)
}

# Filter for INS/DEL only and remove repeat_ids=None
cat("\nFiltering variants (INS/DEL only, removing repeat_ids=None)...\n")
variants_filtered <- variants_df %>%
  filter(SVTYPE %in% c("INS", "DEL")) %>%
  filter(repeat_ids != "None" & !is.na(repeat_ids))

cat("Found", nrow(variants_filtered), "variants with TE annotations\n")

# Function to extract TE families from repeat_ids
extract_te_families <- function(repeat_ids_string) {
  if (is.na(repeat_ids_string) || repeat_ids_string == "None") {
    return(character(0))
  }
  
  annotations <- strsplit(repeat_ids_string, ",")[[1]]
  families <- sapply(annotations, function(x) {
    parts <- strsplit(x, ":")[[1]]
    if (length(parts) >= 1) {
      return(parts[1])
    } else {
      return(NA)
    }
  })
  
  families <- gsub('"', '', families)
  families <- trimws(families)
  return(paste(families, collapse = ";"))
}

# Extract families for each variant
cat("\nExtracting TE families from repeat_ids...\n")
variants_filtered$TE_Families <- sapply(variants_filtered$repeat_ids, extract_te_families)

# Classify each variant into one of the main families
classify_te_family <- function(repeat_ids_string) {
  if (is.na(repeat_ids_string) || repeat_ids_string == "None") {
    return("Other")
  }
  
  line1_patterns <- c("L1", "LINE1", "L1HS", "L1PA", "L1M", "L1ME")
  alu_patterns <- c("Alu", "ALU", "AluJ", "AluS", "AluY")
  sva_patterns <- c("SVA", "SVA_", "SVA-F")
  herv_patterns <- c("HERV", "HER", "LTR", "ERV")
  
  if (any(grepl(paste(line1_patterns, collapse = "|"), repeat_ids_string, ignore.case = TRUE))) {
    return("LINE1")
  } else if (any(grepl(paste(alu_patterns, collapse = "|"), repeat_ids_string, ignore.case = TRUE))) {
    return("Alu")
  } else if (any(grepl(paste(sva_patterns, collapse = "|"), repeat_ids_string, ignore.case = TRUE))) {
    return("SVA")
  } else if (any(grepl(paste(herv_patterns, collapse = "|"), repeat_ids_string, ignore.case = TRUE))) {
    return("HERV")
  }
  
  return("Other")
}

# Apply classification
variants_filtered$Family <- sapply(variants_filtered$repeat_ids, classify_te_family)

# Keep only the main families (remove "Other")
variants_main <- variants_filtered %>%
  filter(Family %in% c("LINE1", "Alu", "SVA", "HERV"))

cat("\nFound", nrow(variants_main), "variants in main TE families\n")

# ============ CALCULATE SUMMARY STATISTICS ============

cat("\nCalculating summary statistics...\n")

# Total TE polymorphisms
total_polymorphisms <- nrow(variants_main)

# Count by family
count_by_family <- table(variants_main$Family)
alu_count <- count_by_family["Alu"]
line1_count <- count_by_family["LINE1"]
sva_count <- count_by_family["SVA"]
herv_count <- count_by_family["HERV"]

# Calculate average per person (mean ALT alleles per individual)
# For each sample, sum the ALT alleles across all variants
alt_alleles_per_sample <- sapply(sample_names, function(sample) {
  sum(variants_main[[paste0(sample, "_ALT")]], na.rm = TRUE)
})

# Calculate mean and standard deviation
mean_per_person <- round(mean(alt_alleles_per_sample), 2)
sd_per_person <- round(sd(alt_alleles_per_sample), 2)
min_per_person <- min(alt_alleles_per_sample)
max_per_person <- max(alt_alleles_per_sample)

# Also calculate by family
alt_by_family <- data.frame(Family = character(), Mean = numeric(), SD = numeric())
for (family in c("LINE1", "Alu", "SVA", "HERV")) {
  family_data <- variants_main %>% filter(Family == family)
  if (nrow(family_data) > 0) {
    family_alt <- sapply(sample_names, function(sample) {
      sum(family_data[[paste0(sample, "_ALT")]], na.rm = TRUE)
    })
    alt_by_family <- rbind(alt_by_family, 
                           data.frame(Family = family, 
                                      Mean = round(mean(family_alt), 2),
                                      SD = round(sd(family_alt), 2)))
  }
}

# ============ CREATE SUMMARY TABLE ============

# Main summary table
summary_table <- data.frame(
  Metric = c(
    "Total TE Polymorphisms",
    "Average per Person (Mean)",
    "Average per Person (SD)",
    "Average per Person (Min)",
    "Average per Person (Max)",
    "",
    "Total Alu",
    "Average Alu per Person",
    "",
    "Total LINE1",
    "Average LINE1 per Person",
    "",
    "Total SVA",
    "Average SVA per Person",
    "",
    "Total HERV",
    "Average HERV per Person"
  ),
  Value = c(
    total_polymorphisms,
    mean_per_person,
    sd_per_person,
    min_per_person,
    max_per_person,
    "",
    alu_count,
    alt_by_family[alt_by_family$Family == "Alu", "Mean"],
    "",
    line1_count,
    alt_by_family[alt_by_family$Family == "LINE1", "Mean"],
    "",
    sva_count,
    alt_by_family[alt_by_family$Family == "SVA", "Mean"],
    "",
    herv_count,
    alt_by_family[alt_by_family$Family == "HERV", "Mean"]
  ),
  stringsAsFactors = FALSE
)

print("\n========================================")
print("SUMMARY TABLE")
print("========================================")
print(summary_table)

# Create a cleaner version with better formatting
summary_clean <- data.frame(
  Category = c("Overall", "Overall", "Overall", "Overall", "Overall",
               "By Family", "By Family", "By Family", "By Family",
               "By Family", "By Family", "By Family", "By Family"),
  Metric = c("Total Polymorphisms", 
             "Mean per Person",
             "SD per Person", 
             "Min per Person",
             "Max per Person",
             "Total Alu", "Mean Alu per Person",
             "Total LINE1", "Mean LINE1 per Person",
             "Total SVA", "Mean SVA per Person",
             "Total HERV", "Mean HERV per Person"),
  Value = c(
    total_polymorphisms,
    mean_per_person,
    sd_per_person,
    min_per_person,
    max_per_person,
    alu_count,
    alt_by_family[alt_by_family$Family == "Alu", "Mean"],
    line1_count,
    alt_by_family[alt_by_family$Family == "LINE1", "Mean"],
    sva_count,
    alt_by_family[alt_by_family$Family == "SVA", "Mean"],
    herv_count,
    alt_by_family[alt_by_family$Family == "HERV", "Mean"]
  )
)

print("\n========================================")
print("CLEAN SUMMARY TABLE")
print("========================================")
print(summary_clean)

# ============ SAVE TABLES ============

# Save as CSV
write.csv(summary_clean, file.path(output_dir, "Summary_Table.csv"), row.names = FALSE)

# Create a formatted text file
sink(file.path(output_dir, "Summary_Table.txt"))

cat("========================================\n")
cat("TRANSPOSON POLYMORPHISM SUMMARY\n")
cat("========================================\n\n")

cat(paste("Analysis Date:", Sys.Date(), "\n"))
cat(paste("Number of samples:", n_samples, "\n"))
cat(paste("Total variants in main families:", total_polymorphisms, "\n\n"))

cat("----------------------------------------\n")
cat("OVERALL STATISTICS\n")
cat("----------------------------------------\n")
cat(paste("Total TE Polymorphisms:", total_polymorphisms, "\n"))
cat(paste("Mean per Person:", mean_per_person, "\n"))
cat(paste("SD per Person:", sd_per_person, "\n"))
cat(paste("Min per Person:", min_per_person, "\n"))
cat(paste("Max per Person:", max_per_person, "\n\n"))

cat("----------------------------------------\n")
cat("BY FAMILY\n")
cat("----------------------------------------\n")
cat(paste("Alu - Total:", alu_count, "\n"))
cat(paste("Alu - Mean per Person:", alt_by_family[alt_by_family$Family == "Alu", "Mean"], "\n\n"))

cat(paste("LINE1 - Total:", line1_count, "\n"))
cat(paste("LINE1 - Mean per Person:", alt_by_family[alt_by_family$Family == "LINE1", "Mean"], "\n\n"))

cat(paste("SVA - Total:", sva_count, "\n"))
cat(paste("SVA - Mean per Person:", alt_by_family[alt_by_family$Family == "SVA", "Mean"], "\n\n"))

cat(paste("HERV - Total:", herv_count, "\n"))
cat(paste("HERV - Mean per Person:", alt_by_family[alt_by_family$Family == "HERV", "Mean"], "\n\n"))

cat("----------------------------------------\n")
cat("DETAILED PER-SAMPLE DATA\n")
cat("----------------------------------------\n")
per_sample_df <- data.frame(
  Sample = sample_names,
  Total_ALT_Alleles = alt_alleles_per_sample
)
print(per_sample_df)

sink()

# Also save per-sample data
write.csv(per_sample_df, file.path(output_dir, "Per_Sample_Data.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("Summary complete!\n")
cat("Files saved to:", output_dir, "\n")
cat("========================================\n")
cat("Files created:\n")
cat("- Summary_Table.csv\n")
cat("- Summary_Table.txt\n")
cat("- Per_Sample_Data.csv\n")
cat("========================================\n")