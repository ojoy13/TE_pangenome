# Load required libraries
library(ggplot2)
library(dplyr)

# Set directories
indir <- "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts_gene_id/perSample"
outdir <- "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts_gene_id/genePlots"

# Create output directory if it doesn't exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Get list of files
file_list <- list.files(indir, pattern = "\\.featureCounts\\.txt$", full.names = TRUE)

# Genes of interest - use the correct versions from your files
genes <- c("ENSG00000075886.11", "ENSG00000198876.13")
gene_names <- c("TUBA3D", "DCAF12")

# Titles for each gene
titles <- c(
  "TUBA3D" = "L1ME3A, AluSc8, AluSx1 deletion in TUBA3D variant Sniffles2.DEL.F4M7C",
  "DCAF12" = "AluY insertion in DCAF12 variant Sniffles2.INS.23MF3"
)

# Genotypes for TUBA3D (from DEL variant F4M7C at chr2:131464150)
genotypes_tuba3d <- c(
  "HG00268" = "0|0",
  "HG03456" = "0|0",
  "HG00358" = "0|0",
  "HG03520" = "0|1",
  "HG01352" = "0|0",
  "HG03807" = "0|0",
  "HG01890" = "1|1",
  "HG04036" = "1|0",
  "HG02059" = "1|0",
  "HG04217" = "1|0",
  "HG02106" = "1|0",
  "NA19129" = "0|1",
  "HG02282" = "0|0",
  "NA19434" = "0|0",
  "HG02769" = "0|0",
  "NA19705" = "1|0",
  "HG02818" = "0|0",
  "NA19836" = "1|0",
  "HG03452" = "0|0",
  "NA20355" = "0|0"
)

# Genotypes for DCAF12 (from INS variant 23MF3 at chr9:34115382)
genotypes_dcaf12 <- c(
  "HG00268" = "1|0",
  "HG03456" = "1|1",
  "HG00358" = "1|1",
  "HG03520" = "1|0",
  "HG01352" = "1|0",
  "HG03807" = "1|1",
  "HG01890" = "1|1",
  "HG04036" = "1|0",
  "HG02059" = "0|0",
  "HG04217" = "1|0",
  "HG02106" = "1|0",
  "NA19129" = "1|0",
  "HG02282" = "1|0",
  "NA19434" = "1|0",
  "HG02769" = "0|1",
  "NA19705" = "1|0",
  "HG02818" = "0|1",
  "NA19836" = "0|0",
  "HG03452" = "0|1",
  "NA20355" = "0|1"
)

# Read all files and extract counts for the two genes
all_counts <- data.frame()

for (file in file_list) {
  # Read file
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get sample name
  sample <- gsub("\\.uniq_2\\.featureCounts\\.txt$", "", basename(file))
  
  # The last column contains the counts
  count_col <- ncol(data)
  
  # Extract counts for the two genes
  for (i in 1:length(genes)) {
    gene_id <- genes[i]
    gene_name <- gene_names[i]
    
    # Check if gene exists in this file
    if (gene_id %in% data$Geneid) {
      count <- data[data$Geneid == gene_id, count_col]
      
      # Assign genotype based on gene
      if (gene_name == "TUBA3D") {
        genotype <- genotypes_tuba3d[sample]
      } else if (gene_name == "DCAF12") {
        genotype <- genotypes_dcaf12[sample]
      }
      
      all_counts <- rbind(all_counts, data.frame(
        Sample = sample,
        Gene = gene_name,
        Genotype = genotype,
        Counts = as.numeric(count)
      ))
    }
  }
}

# Check what we collected
print("All counts collected:")
print(head(all_counts))
print(str(all_counts))

# Create separate violin plots for each gene
for (gene in gene_names) {
  # Subset data for this gene
  gene_data <- all_counts[all_counts$Gene == gene, ]
  
  # Check if we have data
  print(paste("Data for", gene, ":", nrow(gene_data), "rows"))
  print(head(gene_data))
  
  # Remove any rows with NA genotype or NA counts
  gene_data <- gene_data[!is.na(gene_data$Genotype) & !is.na(gene_data$Counts), ]
  
  # Order genotypes properly
  gene_data$Genotype <- factor(gene_data$Genotype, levels = c("0|0", "0|1", "1|0", "1|1"))
  
  # Skip if no data
  if (nrow(gene_data) == 0) {
    print(paste("No data for", gene, "- skipping"))
    next
  }
  
  # Calculate summary statistics for each genotype
  summary_stats <- gene_data %>%
    group_by(Genotype) %>%
    summarise(
      Mean = mean(Counts, na.rm = TRUE),
      SD = sd(Counts, na.rm = TRUE),
      N = n()
    )
  
  # Print summary
  cat("\n", gene, "summary by genotype:\n")
  print(summary_stats)
  cat("\n")
  
  # Create violin plot with boxplot overlay
  p <- ggplot(gene_data, aes(x = Genotype, y = Counts)) +
    # Violin plot
    geom_violin(fill = "white", color = "black", width = 0.8) +
    # Boxplot overlay
    geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
    # Jitter points (optional)
    geom_jitter(width = 0.1, alpha = 0.3, size = 1, color = "black") +
    # Labels and theme
    labs(
      title = titles[gene],
      x = "Genotype",
      y = "Feature Counts"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.line = element_line(color = "black")
    )
  
  # Save plots
  ggsave(file.path(outdir, paste0(gene, "_violin_with_boxplot.pdf")), p, width = 8, height = 6)
  ggsave(file.path(outdir, paste0(gene, "_violin_with_boxplot.png")), p, width = 8, height = 6, dpi = 300)
  
  # Also save the summary statistics
  write.csv(summary_stats, file.path(outdir, paste0(gene, "_summary_stats.csv")), row.names = FALSE)
}

print("Done! Check output directory for plots and summary statistics.")