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

# Initialize Location column
variants_main$Location <- "intergenic"
variants_main$Location[exon_hits] <- "exonic"

# Variants in genes but not exons = intronic
intronic_hits <- setdiff(gene_hits, exon_hits)
variants_main$Location[intronic_hits] <- "intronic"

# ============ ENHANCER ANNOTATION ============
cat("\nLoading enhancer annotations...\n")
enhancer_file <- "/scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/enhancers/pang21_dreg_enhancers.bed"

# Read the enhancer BED file
enhancer_df <- read.table(enhancer_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
cat("Total enhancer regions:", nrow(enhancer_df), "\n")

# Let's examine the structure of the enhancer file
cat("\nEnhancer file structure:\n")
print(head(enhancer_df, 3))

# Column 5 appears to be the variant ID from the VCF
# Let's extract variant IDs from the enhancer file
enhancer_variant_ids <- unique(enhancer_df[,5])
cat("\nUnique variant IDs in enhancer file:", length(enhancer_variant_ids), "\n")

# Show some examples
cat("Example enhancer variant IDs:\n")
print(head(enhancer_variant_ids, 10))

# Check if any of these IDs match our variants
matching_ids <- variants_main$ID[variants_main$ID %in% enhancer_variant_ids]
cat("\nMatching variant IDs found:", length(matching_ids), "\n")

# Show some matching examples if any
if (length(matching_ids) > 0) {
  cat("Matching variant IDs:\n")
  print(head(matching_ids, 10))
} else {
  cat("No matching variant IDs found!\n")
  cat("Let's check if there's a formatting issue...\n")
  cat("First few variant IDs from VCF:\n")
  print(head(variants_main$ID, 10))
  cat("\nFirst few variant IDs from enhancers:\n")
  print(head(enhancer_variant_ids, 10))
}

# Alternative: Try matching by position and chromosome instead of ID
# Create GRanges for enhancers (using columns 1-3 for chrom, start, end)
cat("\nTrying to match by genomic coordinates instead...\n")

# The enhancer BED might have: chr, start, end, name, variant_id, ... 
# Let's use columns 1-3 for coordinates
enhancer_gr <- GRanges(
  seqnames = enhancer_df[,1],
  ranges = IRanges(start = enhancer_df[,2], end = enhancer_df[,3])
)

# Find overlaps between variants and enhancers
enhancer_overlaps <- findOverlaps(variant_gr, enhancer_gr)
enhancer_hits <- unique(queryHits(enhancer_overlaps))

cat("Variants overlapping enhancer regions by coordinates:", length(enhancer_hits), "\n")

if (length(enhancer_hits) > 0) {
  # Update Location for enhancer hits
  variants_main$Location[enhancer_hits] <- "enhancer"
  cat("Updated", length(enhancer_hits), "variants to 'enhancer'\n")
} else {
  cat("No variants overlap enhancer regions by coordinates.\n")
  cat("Let's check coordinate ranges...\n")
  cat("Variant coordinate range:\n")
  print(range(variants_main$POS))
  cat("\nEnhancer coordinate range:\n")
  print(range(enhancer_df[,2:3]))
}

# Also check if enhancers are on the same chromosomes
cat("\nChromosomes in variants:\n")
print(unique(seqnames(variant_gr))[1:10])
cat("\nChromosomes in enhancers:\n")
print(unique(enhancer_df[,1])[1:10])

# If still no hits, maybe the enhancer file needs 0-based to 1-based conversion
if (length(enhancer_hits) == 0) {
  cat("\nTrying with 0-based to 1-based conversion for enhancers...\n")
  enhancer_gr_0based <- GRanges(
    seqnames = enhancer_df[,1],
    ranges = IRanges(start = enhancer_df[,2] + 1, end = enhancer_df[,3])
  )
  enhancer_overlaps_0based <- findOverlaps(variant_gr, enhancer_gr_0based)
  enhancer_hits_0based <- unique(queryHits(enhancer_overlaps_0based))
  cat("Variants overlapping enhancers with 0-based conversion:", length(enhancer_hits_0based), "\n")
  
  if (length(enhancer_hits_0based) > 0) {
    variants_main$Location[enhancer_hits_0based] <- "enhancer"
    cat("Updated", length(enhancer_hits_0based), "variants to 'enhancer'\n")
  }
}

cat("\nFinal Location breakdown:\n")
print(table(variants_main$Location))
cat("\nLocation breakdown by family:\n")
print(table(variants_main$Family, variants_main$Location))

# ============ CREATE PLOT ============

# Colorblind-friendly palette for locations (including enhancer)
location_colors <- c(
  "intergenic" = "#999999",  # Gray
  "intronic" = "#66CCEE",    # Light blue
  "exonic" = "#EE6677",      # Red/pink
  "enhancer" = "#44BB99"     # Teal/green
)

# Set the order of families (Alu, LINE1, SVA, HERV)
family_order <- c("Alu", "LINE1", "SVA", "HERV")

# Colorblind-friendly palette for TE families (in desired order)
family_colors <- c(
  "Alu" = "#56B4E9",    # Sky blue
  "LINE1" = "#E69F00",  # Orange
  "SVA" = "#009E73",    # Bluish green
  "HERV" = "#CC79A7"    # Magenta
)

# Prepare data for plotting - grouped bars (dodged) with all 4 locations
plot_data <- variants_main %>%
  group_by(Location, Family) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  complete(Location = c("intergenic", "intronic", "exonic", "enhancer"), 
           Family = family_order, 
           fill = list(Count = 0))

# Order locations for consistent display
plot_data$Location <- factor(plot_data$Location, 
                             levels = c("intergenic", "intronic", "exonic", "enhancer"))

# Order Family according to desired order
plot_data$Family <- factor(plot_data$Family, levels = family_order)

# Calculate max for y-axis
max_count <- max(plot_data$Count)
y_lim <- max_count * 1.15

# Create publication-ready grouped bar plot
cat("\nCreating location plot with enhancers...\n")

plot_location <- ggplot(plot_data, aes(x = Location, y = Count, fill = Family)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = family_colors) +
  scale_x_discrete(labels = c("intergenic" = "Intergenic", 
                              "intronic" = "Intronic", 
                              "exonic" = "Exonic",
                              "enhancer" = "Enhancer")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     limits = c(0, y_lim)) +
  labs(
    title = "Genomic Location of Transposons",
    subtitle = "By TE Family",
    x = "Genomic Location",
    y = "Count",
    fill = "Transposon Family"
  ) +
  theme_classic() +
  theme(
    # Text sizes for publication
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    # Legend inside plot
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
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
ggsave(file.path(output_dir, "transposon_location_by_family_with_enhancers.png"), 
       plot_location, width = 14, height = 8, dpi = 600)

ggsave(file.path(output_dir, "transposon_location_by_family_with_enhancers.pdf"), 
       plot_location, width = 14, height = 8, device = cairo_pdf)

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
                              "exonic" = "Exonic",
                              "enhancer" = "Enhancer")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     limits = c(0, y_lim_pct)) +
  labs(
    title = "Genomic Location of Transposon",
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
ggsave(file.path(output_dir, "transposon_location_by_family_percentage_with_enhancers.png"), 
       plot_location_pct, width = 14, height = 8, dpi = 600)

ggsave(file.path(output_dir, "transposon_location_by_family_percentage_with_enhancers.pdf"), 
       plot_location_pct, width = 14, height = 8, device = cairo_pdf)

# Create summary table
location_summary <- variants_main %>%
  group_by(Location, Family) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  complete(Location = c("intergenic", "intronic", "exonic", "enhancer"), 
           Family = family_order, 
           fill = list(Count = 0)) %>%
  pivot_wider(names_from = Family, values_from = Count, values_fill = 0) %>%
  mutate(Total = Alu + LINE1 + SVA + HERV)

# Reorder locations
location_summary$Location <- factor(location_summary$Location, 
                                    levels = c("intergenic", "intronic", "exonic", "enhancer"))
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
write.csv(location_summary, file.path(output_dir, "location_by_family_summary_with_enhancers.csv"), row.names = FALSE)
write.csv(location_summary_pct, file.path(output_dir, "location_by_family_summary_with_enhancers_pct.csv"), row.names = FALSE)

# Save annotated data
variants_main_export <- variants_main %>%
  select(-TE_Families)
write.csv(variants_main_export, file.path(output_dir, "annotated_variants_with_enhancers.csv"), row.names = FALSE)

# Create README
readme_content <- c(
  "# Transposon Genomic Location Analysis (with Enhancers)",
  "",
  paste0("Generated on: ", Sys.Date()),
  "",
  "## Files in this directory:",
  "",
  "1. **transposon_location_by_family_with_enhancers.png** - Grouped bar plot with enhancers (PNG, 600 DPI)",
  "2. **transposon_location_by_family_with_enhancers.pdf** - Grouped bar plot with enhancers (Vector PDF)",
  "3. **transposon_location_by_family_percentage_with_enhancers.png** - Percentage version (PNG, 600 DPI)",
  "4. **transposon_location_by_family_percentage_with_enhancers.pdf** - Percentage version (Vector PDF)",
  "5. **location_by_family_summary_with_enhancers.csv** - Summary table with counts",
  "6. **location_by_family_summary_with_enhancers_pct.csv** - Summary table with percentages",
  "7. **annotated_variants_with_enhancers.csv** - All variants with location annotation",
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
  "- Enhancer: Overlapping an enhancer region (may also overlap other categories)",
  "",
  "## Counts by family:",
  paste0("- Alu: ", sum(variants_main$Family == "Alu")),
  paste0("- LINE1: ", sum(variants_main$Family == "LINE1")),
  paste0("- SVA: ", sum(variants_main$Family == "SVA")),
  paste0("- HERV: ", sum(variants_main$Family == "HERV")),
  "",
  "## Counts by location:",
  paste0("- Intergenic: ", sum(variants_main$Location == "intergenic")),
  paste0("- Intronic: ", sum(variants_main$Location == "intronic")),
  paste0("- Exonic: ", sum(variants_main$Location == "exonic")),
  paste0("- Enhancer: ", sum(variants_main$Location == "enhancer"))
)

writeLines(readme_content, file.path(output_dir, "README_with_enhancers.txt"))

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("All outputs saved to:", output_dir, "\n")
cat("========================================\n")
cat("Files created:\n")
cat(list.files(output_dir), sep = "\n")
cat("\n========================================\n")
cat("Location summary:\n")
print(location_summary)
cat("========================================\n")