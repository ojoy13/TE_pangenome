# Load required libraries
library(tidyverse)
library(vcfR)
library(ggplot2)
library(scales)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)

# Create output directory
output_dir <- "TE_family_location"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created directory:", output_dir, "\n")
}

# Read the VCF file
vcf_file <- "/scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/pangenome_21_dedup.vcf"
vcf <- read.vcfR(vcf_file)

# Extract INFO fields
info_df <- as.data.frame(vcf@fix[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")])

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

# Extract SVTYPE, repeat_ids, and END
cat("Parsing VCF INFO fields...\n")
variants_df <- data.frame(
  CHROM = info_df$CHROM,
  POS = as.numeric(info_df$POS),
  ID = info_df$ID,
  SVTYPE = sapply(info_df$INFO, function(x) extract_info_value(x, "SVTYPE")),
  END = as.numeric(sapply(info_df$INFO, function(x) extract_info_value(x, "END"))),
  repeat_ids = sapply(info_df$INFO, function(x) extract_info_value(x, "repeat_ids")),
  stringsAsFactors = FALSE
)

# For variants without END, use POS (for INS) or POS+1 (for DEL)
variants_df$END <- ifelse(is.na(variants_df$END), 
                          ifelse(variants_df$SVTYPE == "INS", variants_df$POS, variants_df$POS + 1),
                          variants_df$END)

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

# ============ GENOMIC LOCATION ANNOTATION ============
cat("\nLoading GTF annotation file...\n")
gtf_file <- "/scratch/Users/olde5615/data/gencode.v38.ComprehensiveAnnotation.Main.gtf"
gtf <- import.gff(gtf_file, format = "gtf")

# Keep only gene and transcript annotations
gtf_genes <- gtf[gtf$type == "gene"]
gtf_exons <- gtf[gtf$type == "exon"]

# Create GRanges for variants
cat("\nCreating variant GRanges...\n")
variant_gr <- GRanges(
  seqnames = variants_main$CHROM,
  ranges = IRanges(start = variants_main$POS, end = variants_main$END),
  mcols = DataFrame(
    ID = variants_main$ID,
    SVTYPE = variants_main$SVTYPE,
    Family = variants_main$Family,
    repeat_ids = variants_main$repeat_ids
  )
)

# Remove any variants with invalid chromosomes (e.g., chrM, chrUn)
valid_chroms <- paste0("chr", c(1:22, "X", "Y"))
variant_gr <- variant_gr[seqnames(variant_gr) %in% valid_chroms]
cat("Using", length(variant_gr), "variants on valid chromosomes\n")

# Annotate genomic locations
cat("\nAnnotating genomic locations...\n")

# Create exonic regions
exon_gr <- GRanges(
  seqnames = gtf_exons@seqnames,
  ranges = gtf_exons@ranges,
  mcols = DataFrame(gene_id = gtf_exons$gene_id)
)

# Find overlaps with exons
exon_overlaps <- findOverlaps(variant_gr, exon_gr)
exon_hits <- unique(queryHits(exon_overlaps))

# Create intronic regions (gene bodies minus exons)
gene_gr <- GRanges(
  seqnames = gtf_genes@seqnames,
  ranges = gtf_genes@ranges,
  mcols = DataFrame(gene_id = gtf_genes$gene_id)
)

# Find variants in genes (including introns)
gene_overlaps <- findOverlaps(variant_gr, gene_gr)
gene_hits <- unique(queryHits(gene_overlaps))

# Variants in exons
variants_main$Location <- "intergenic"
variants_main$Location[exon_hits] <- "exonic"

# Variants in genes but not exons = intronic
intronic_hits <- setdiff(gene_hits, exon_hits)
variants_main$Location[intronic_hits] <- "intronic"

cat("\nLocation breakdown:\n")
print(table(variants_main$Location))
cat("\nLocation breakdown by family:\n")
print(table(variants_main$Family, variants_main$Location))

# ============ CREATE PLOT ============

# Colorblind-friendly palette for TE families (in desired order)
family_colors <- c(
  "Alu" = "#56B4E9",    # Sky blue
  "LINE1" = "#E69F00",  # Orange
  "SVA" = "#009E73",    # Bluish green
  "HERV" = "#CC79A7"    # Magenta
)

# Set the order of families (Alu, LINE1, SVA, HERV)
family_order <- c("Alu", "LINE1", "SVA", "HERV")

# Prepare data for plotting - grouped bars (dodged)
plot_data <- variants_main %>%
  group_by(Location, Family) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  complete(Location = c("intergenic", "intronic", "exonic"), 
           Family = family_order, 
           fill = list(Count = 0))

# Order locations for consistent display
plot_data$Location <- factor(plot_data$Location, 
                             levels = c("intergenic", "intronic", "exonic"))

# Order Family according to desired order
plot_data$Family <- factor(plot_data$Family, levels = family_order)

# Calculate max for y-axis
max_count <- max(plot_data$Count)
y_lim <- max_count * 1.15

# Create publication-ready grouped bar plot
cat("\nCreating location plot...\n")

plot_location <- ggplot(plot_data, aes(x = Location, y = Count, fill = Family)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = family_colors) +
  scale_x_discrete(labels = c("intergenic" = "Intergenic", 
                              "intronic" = "Intronic", 
                              "exonic" = "Exonic")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     limits = c(0, y_lim)) +
  labs(
    title = "Genomic Location of TE Polymorphisms",
    subtitle = "By TE Family",
    x = "Genomic Location",
    y = "Count",
    fill = "Transposon Family"
  ) +
  theme_classic() +
  theme(
    # Text sizes for publication
    plot.title = element_text(size = 25, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    # Legend inside plot
    legend.position = c(0.85, 0.75),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    legend.spacing.y = unit(0.3, "cm"),
    legend.key.size = unit(0.8, "cm"),
    legend.key.spacing = unit(0.3, "cm"),
    legend.box = "vertical",
    # Plot margins
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
  ) +
  # Add count labels on bars
  geom_text(aes(label = ifelse(Count > 0, Count, ""), 
                group = Family),
            position = position_dodge(width = 0.8), 
            vjust = -0.5, 
            size = 6,
            fontface = "bold")

print(plot_location)

# Save plots
ggsave(file.path(output_dir, "transposon_location_by_family.png"), 
       plot_location, width = 10.5, height = 8, dpi = 600)

ggsave(file.path(output_dir, "transposon_location_by_family.pdf"), 
       plot_location, width = 10.5, height = 8, device = cairo_pdf)

# Create a version with percentage instead of counts
plot_data_pct <- plot_data %>%
  group_by(Location) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Calculate max percentage for y-axis
max_pct <- max(plot_data_pct$Percentage)
y_lim_pct <- max_pct * 1.15

plot_location_pct <- ggplot(plot_data_pct, aes(x = Location, y = Percentage, fill = Family)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = family_colors) +
  scale_x_discrete(labels = c("intergenic" = "Intergenic", 
                              "intronic" = "Intronic", 
                              "exonic" = "Exonic")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     limits = c(0, y_lim_pct)) +
  labs(
    title = "Genomic Location of Transposon Insertions and Deletions",
    subtitle = "By TE Family (Percentage)",
    x = "Genomic Location",
    y = "Percentage (%)",
    fill = "Transposon Family"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    legend.spacing.y = unit(0.3, "cm"),
    legend.key.size = unit(0.8, "cm"),
    legend.key.spacing = unit(0.3, "cm"),
    legend.box = "vertical",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
  ) +
  geom_text(aes(label = ifelse(Percentage > 0, paste0(round(Percentage, 1), "%"), ""), 
                group = Family),
            position = position_dodge(width = 0.8), 
            vjust = -0.5, 
            size = 6,
            fontface = "bold")

print(plot_location_pct)

# Save percentage version
ggsave(file.path(output_dir, "transposon_location_by_family_percentage.png"), 
       plot_location_pct, width = 12, height = 8, dpi = 600)

ggsave(file.path(output_dir, "transposon_location_by_family_percentage.pdf"), 
       plot_location_pct, width = 12, height = 8, device = cairo_pdf)

# Create summary table
location_summary <- variants_main %>%
  group_by(Location, Family) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  complete(Location = c("intergenic", "intronic", "exonic"), 
           Family = family_order, 
           fill = list(Count = 0)) %>%
  pivot_wider(names_from = Family, values_from = Count, values_fill = 0) %>%
  mutate(Total = Alu + LINE1 + SVA + HERV)

# Reorder locations
location_summary$Location <- factor(location_summary$Location, 
                                    levels = c("intergenic", "intronic", "exonic"))
location_summary <- location_summary %>% arrange(Location)

print("\nLocation summary by family:")
print(location_summary)

# Also create a table with percentages
location_summary_pct <- location_summary %>%
  mutate(
    Alu_pct = round(Alu / Total * 100, 1),
    LINE1_pct = round(LINE1 / Total * 100, 1),
    SVA_pct = round(SVA / Total * 100, 1),
    HERV_pct = round(HERV / Total * 100, 1)
  )

print("\nLocation summary with percentages:")
print(location_summary_pct)

# Save summary tables
write.csv(location_summary, file.path(output_dir, "location_by_family_summary.csv"), row.names = FALSE)
write.csv(location_summary_pct, file.path(output_dir, "location_by_family_summary_with_pct.csv"), row.names = FALSE)

# Save annotated data
variants_main_export <- variants_main %>%
  select(-TE_Families)
write.csv(variants_main_export, file.path(output_dir, "annotated_variants_with_location.csv"), row.names = FALSE)

# Create README
readme_content <- c(
  "# Transposon Genomic Location Analysis",
  "",
  paste0("Generated on: ", Sys.Date()),
  "",
  "## Files in this directory:",
  "",
  "1. **transposon_location_by_family.png** - Grouped bar plot (PNG, 600 DPI)",
  "2. **transposon_location_by_family.pdf** - Grouped bar plot (Vector PDF)",
  "3. **transposon_location_by_family_percentage.png** - Percentage version (PNG, 600 DPI)",
  "4. **transposon_location_by_family_percentage.pdf** - Percentage version (Vector PDF)",
  "5. **location_by_family_summary.csv** - Summary table with counts",
  "6. **location_by_family_summary_with_pct.csv** - Summary table with percentages",
  "7. **annotated_variants_with_location.csv** - All variants with location annotation",
  "",
  "## Color palette (colorblind-friendly):",
  "- Alu: Sky Blue (#56B4E9)",
  "- LINE1: Orange (#E69F00)",
  "- SVA: Bluish Green (#009E73)",
  "- HERV: Magenta (#CC79A7)",
  "",
  "## Filtering parameters:",
  "- SVTYPE: INS or DEL only",
  "- Removed variants with repeat_ids=None",
  "- Main TE families only (LINE1, Alu, SVA, HERV)",
  "",
  paste0("Total variants: ", nrow(variants_main)),
  "",
  "## Location definitions:",
  "- Intergenic: Outside of gene boundaries",
  "- Intronic: Within a gene but not in an exon",
  "- Exonic: Overlapping an exon",
  "",
  "## Counts by family:",
  paste0("- Alu: ", sum(variants_main$Family == "Alu")),
  paste0("- LINE1: ", sum(variants_main$Family == "LINE1")),
  paste0("- SVA: ", sum(variants_main$Family == "SVA")),
  paste0("- HERV: ", sum(variants_main$Family == "HERV"))
)

writeLines(readme_content, file.path(output_dir, "README_location.txt"))

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("All outputs saved to:", output_dir, "\n")
cat("========================================\n")
cat("Files created:\n")
cat(list.files(output_dir), sep = "\n")
cat("\n========================================\n")
cat("Location summary:\n")
print(location_summary)
cat("\n========================================\n")