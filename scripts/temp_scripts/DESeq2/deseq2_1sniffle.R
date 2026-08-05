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

# Format metadata
meta <- read.csv("/scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/21_sample_metadata.csv")
sniffle <- read.csv("/scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/deseq2_metadata.csv")

# Filter for your favorite Sniffle (you can change this ID)
onesniffle <- sniffle %>% filter(sniffles_id == "Sniffles2.INS.262M1D")
onesniffle <- onesniffle %>% select(-chromosome, -position, -sniffles_id, -te_types)
onesniffleT <- t(onesniffle)
colnames(onesniffleT) <- c("TE_group")
onesniffleT <- as.data.frame(onesniffleT)  
onesniffleT$bamFilename <- rownames(onesniffleT)
onesniffleT <- onesniffleT %>% separate_wider_delim(cols = bamFilename, delim = ".", names = c("SampleID", "merge", "bam")) 
onesniffleT <- onesniffleT %>% select(TE_group, SampleID)

meta2 <- merge(meta, onesniffleT, by = "SampleID")

# Filter metadata to only samples that exist in counts
meta3 <- as.data.frame(meta2) %>% filter(SampleID %in% colnames(counts_df))

# Prepare count matrix for DESeq2
rownames(counts_df) <- counts_df$Geneid
counts_df <- counts_df %>% select(-Geneid)
counts_df <- counts_df %>% select(all_of(meta3$SampleID))  # Reorder columns

# Set TE_group as factor
meta3$TE_group <- as.factor(meta3$TE_group)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_df, 
                              colData = meta3, 
                              design = ~ TE_group)

# Run DESeq2
DEdds <- DESeq(dds)

# Check size factors
sizeFactors(DEdds)

# Get results:
res <- results(DEdds)
head(res)
# Visualization: Before and after normalization
# Raw counts
outdir <- "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts_gene_id/DESeq2_results/"

# Make sure directory exists
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 1. Raw counts violin plot
allcounts <- as.data.frame(counts(DEdds))
allcountslong <- allcounts %>% gather(key = "sample", value = "signal")

p1 <- ggplot(allcountslong, aes(x = sample, y = signal)) + 
  geom_violin(trim = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) + 
  scale_y_continuous(trans = 'log2') +
  ggtitle("Raw Counts Distribution")
ggsave(paste0(outdir, "RNA_raw_counts_violin_linear_gene.pdf"), p1, width = 14, height = 6)

# 2. Normalized counts violin plot
normcounts <- as.data.frame(counts(DEdds, normalize = TRUE))
normcountslong <- normcounts %>% gather(key = "sample", value = "signal")

p2 <- ggplot(normcountslong, aes(x = sample, y = signal)) + 
  geom_violin(trim = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) + 
  scale_y_continuous(trans = 'log2') +
  ggtitle("Normalized Counts Distribution")
ggsave(paste0(outdir, "RNA_normalized_counts_violin_linear_gene.pdf"), p2, width = 14, height = 6)

# 3. Dispersion plot
pdf(paste0(outdir, "RNA_dispersion_estimates_linear_gene.pdf"), width = 8, height = 6)
plotDispEsts(DEdds)
title("Dispersion Estimates")
dev.off()

# 4. PCA plot
vsd <- vst(DEdds, blind = FALSE)
p3 <- plotPCA(vsd, intgroup = "TE_group") +
  theme_minimal() +
  ggtitle("PCA Plot - RNA Samples")
ggsave(paste0(outdir, "RNA_PCA_plot_linear_gene.pdf"), p3, width = 8, height = 6)

# 5. MA Plot
pdf(paste0(outdir, "RNA_MA_plot_linear_gene.pdf"), width = 8, height = 6)
plotMA(res)
title("MA Plot")
dev.off()

# 6. Save results summary to text file
sink(paste0(outdir, "RNA_DESeq2_results_summary_linear_gene.txt"))
print("Results summary:")
summary(res)
print("\nNumber of significant genes (padj < 0.05):")
sum(res$padj < 0.05, na.rm = TRUE)
print("\nNumber of significant genes (padj < 0.01):")
sum(res$padj < 0.01, na.rm = TRUE)
print("\nTop 20 genes by padj:")
head(res[order(res$padj), ], 20)
sink()

# 7. Save significant genes to CSV
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
sig_genes <- res_df %>% 
  filter(padj < 0.05, !is.na(padj)) %>%
  arrange(padj)
write.csv(sig_genes, paste0(outdir, "RNA_significant_geneIDs_padj0.05_linear_gene.csv"), row.names = FALSE)

write.csv(res_df, paste0(outdir, "RNA_all_geneIDs_linear_gene.csv"), row.names = FALSE)

# 8. Sample size information
sample_counts <- table(meta3$TE_group)
write.csv(as.data.frame(sample_counts), paste0(outdir, "RNA_sample_sizes_linear_gene.csv"))

# Plot dispersion estimates
plotDispEsts(DEdds)

# Get results - you need to specify which levels to compare
# First check what levels exist in your TE_group
levels(meta3$TE_group)

# Set your comparison groups (replace with your actual group names)
# Example: if TE_group has levels "0" and "1"
sample1 <- "0"  # Control/Reference group
sample2 <- "1"  # Test group

# OR if you have specific group names like "D21_thirtyseven" etc.
# sample1 <- "D21_thirtyseven"
# sample2 <- "D21_fourtytwo"

# Get results
res <- results(DEdds, contrast = c("TE_group", sample1, sample2))

# Save results
fileroot <- "all_samples_Sniffles2.INS.262M1D_linear_gene"
write.csv(res, paste0(outdir, fileroot, "_", sample1, "_vs_", sample2, ".results.csv"))

# Look at results
summary(res)
plotMA(res)


# Load GTF for gene annotations
gtf <- "/Shares/CL_Shared/db/genomes/hg38/annotations/gencode.v38.ComprehensiveAnnotation.Main.gtf"

# Read the GTF file
gtf_df <- read.table(gtf, 
                            sep = "\t", 
                            header = FALSE, 
                            quote = "",
                            stringsAsFactors = FALSE,
                            comment.char = "#")

# Add column names for GTF (standard GTF format)
colnames(gtf_df) <- c("seqname", "source", "feature", "start", "end", 
                             "score", "strand", "frame", "attribute")
gtf_genes <- gtf_df %>% filter(feature == "transcript")

# Read all transcripts
# Base R solution
extract_field <- function(attr, field) {
  pattern <- paste0(field, " ([^;]+)")
  match <- regexpr(pattern, attr)
  ifelse(match > 0, 
         substr(attr, match + nchar(field) + 1, 
                match + attr(match, "match.length") - 1),
         NA)
}

gtf_transcripts$transcript_id <- extract_field(gtf_transcripts$attribute, "gene_id")
gtf_transcripts$gene_name <- extract_field(gtf_transcripts$attribute, "gene_name")

# Check
head(gtf_transcripts[, c("gene_id", "gene_name")])

# Get unique
transcript_unique <- unique(gtf_transcripts[!is.na(gtf_genes$gene_id), 
                                            c("transcript_id", "gene_name")])
nrow(transcript_unique)

# Read "main" transcripts
# Remove version numbers for matching
transcript_unique$transcript_base <- gsub("\\..*$", "", transcript_unique$transcript_id)

# Get your count transcript IDs (without versions)
counts_base <- gsub("\\..*$", "", rownames(counts_df))

# Find transcripts that exist in your count data
main_base <- intersect(transcript_unique$transcript_base, counts_base)
length(main_base)  # This should be much smaller (~20-30k)

# Create mapping with full IDs
main_mapping <- transcript_unique[transcript_unique$transcript_base %in% main_base, 
                                  c("transcript_base", "transcript_id", "gene_name")]
main_mapping <- unique(main_mapping)

# Add the full transcript IDs from your counts (with version numbers)
counts_df_temp <- data.frame(
  transcript_base = counts_base,
  full_transcript_id = rownames(counts_df),
  stringsAsFactors = FALSE
)
counts_df_temp <- unique(counts_df_temp)

main_mapping <- merge(main_mapping, counts_df_temp, by = "transcript_base")
nrow(main_mapping)  # Your number of main transcripts

# Filter your results
res_df <- as.data.frame(res)
res_df$transcript_id <- rownames(res_df)

# Keep only main transcripts
res_main <- res_df[res_df$transcript_id %in% main_mapping$full_transcript_id, ]

# Add gene names
res_main <- merge(res_main, main_mapping[, c("full_transcript_id", "gene_name", "transcript_id")], 
                  by.x = "transcript_id", by.y = "full_transcript_id", 
                  all.x = TRUE)

# Sort by significance
res_main <- res_main[order(res_main$padj), ]

# Results summary
cat("\n=== Main Transcripts Analysis ===\n")
cat("Total transcripts in Gencode:", nrow(transcript_unique), "\n")
cat("Transcripts in your count data:", length(unique(counts_base)), "\n")
cat("Main transcripts analyzed:", nrow(res_main), "\n")
cat("Significant main transcripts (padj < 0.05):", 
    sum(res_main$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant main transcripts (padj < 0.01):", 
    sum(res_main$padj < 0.01, na.rm = TRUE), "\n")

# View top significant main transcripts
cat("\n=== Top 20 Significant Main Transcripts ===\n")
print(head(res_main[, c("gene_name", "transcript_id", "log2FoldChange", "padj")], 20))

# Save results
write.csv(res_main, paste0(outdir, "RNA_main_transcripts_all_results.csv"), 
          row.names = FALSE)

# Save only significant ones
sig_main <- res_main[!is.na(res_main$padj) & res_main$padj < 0.05, ]
sig_main <- sig_main[order(sig_main$padj), ]
print(head(sig_main[, c("gene_name", "transcript_id", "log2FoldChange", "padj")], 20))
write.csv(sig_main, paste0(outdir, "RNA_main_transcripts_significant.csv"), 
          row.names = FALSE)

# Create a summary plot
library(ggplot2)

# Histogram of padj for main transcripts
ggplot(res_main, aes(x = -log10(padj))) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed") +
  labs(title = "Distribution of Significance for Main Transcripts",
       x = "-log10(Adjusted P-value)",
       y = "Count") +
  theme_bw()

ggsave(paste0(outdir, "RNA_main_transcripts_significance_distribution.pdf"), 
       width = 8, height = 6)


# re-run deseq2 without chrY
# First, add chromosome info to gtf_transcripts
gtf_transcripts$chr <- gtf_df$seqname[match(rownames(gtf_transcripts), rownames(gtf_df))]

# Get unique Y-chromosome genes
chrY_genes <- unique(gtf_transcripts$gene_name[gtf_transcripts$chr == "chrY" & 
                                                 !is.na(gtf_transcripts$gene_name)])
length(chrY_genes)

# Combine with your results-based list
y_genes_from_results <- c("UTY", "TXLNGY", "KDM5D", "USP9Y", "DDX3Y", "TTTY10")
y_genes <- c("UTY", "TXLNGY", "KDM5D", "USP9Y", "DDX3Y", "TTTY10", 
             "RPS4Y1", "RPS4Y2", "EIF1AY", "ZFY", "CYorf15A", "CYorf15B",
             "TMSB4Y", "NLGN4Y", "PCDH11Y", "TGIF2LY", "TSPY1", "TSPY2",
             "TSPY3", "TSPY4", "TSPY8", "TSPY9P", "VCY", "XKRY", "HSFY1",
             "HSFY2", "PRKY", "RBMY1A1", "RBMY1B", "RBMY1C", "RBMY1D", 
             "RBMY1E", "RBMY1F", "RBMY1J", "BPY2", "CDY1", "CDY2", 
             "DAZ1", "DAZ2", "DAZ3", "DAZ4")
all_y_genes <- unique(c(y_genes, y_genes_from_results, chrY_genes))
length(all_y_genes)

# Now filter again
res_filtered <- res_main[!res_main$gene_name %in% all_y_genes, ]
res_filtered <- res_filtered[order(res_filtered$padj), ]

# Step 6: View top results without Y genes
cat("\n=== Top 20 Results After Removing Y-chromosome Genes ===\n")
print(head(res_filtered[, c("gene_name", "log2FoldChange", "pvalue", "padj")], 20))

# Step 7: Check for significant genes
sig_filtered <- res_filtered[!is.na(res_filtered$padj) & res_filtered$padj < 0.05, ]
cat("\nSignificant non-Y genes (padj < 0.05):", nrow(sig_filtered), "\n")

if(nrow(sig_filtered) > 0) {
  print(sig_filtered[, c("gene_name", "log2FoldChange", "padj")])
} else {
  cat("\nNo significant genes at padj < 0.05\n")
  cat("Top 20 by raw p-value:\n")
  print(head(res_filtered[!is.na(res_filtered$pvalue), 
                          c("gene_name", "log2FoldChange", "pvalue")], 20))
}

# Step 8: Save filtered results
write.csv(res_filtered, 
          paste0(outdir, "RNA_results_without_Y_genes.csv"), 
          row.names = FALSE)

# Step 9: Create volcano plot without Y genes
library(ggplot2)

ggplot(res_filtered, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.6, size = 1, color = "steelblue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "DESeq2 Results - Y Chromosome Genes Removed",
       x = "Log2 Fold Change", 
       y = "-Log10 P-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(outdir, "RNA_volcano_plot_no_Y.pdf"), width = 8, height = 6)

# Step 10: Create a bar plot of top 20 genes by p-value
top20 <- head(res_filtered[!is.na(res_filtered$pvalue), ], 20)

ggplot(top20, aes(x = reorder(gene_name, -log10(pvalue)), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = "Top 20 Genes by P-value (Y-chromosome Removed)",
       x = "Gene Name",
       y = "-log10(P-value)") +
  theme_bw()

ggsave(paste0(outdir, "RNA_top20_genes_no_Y.pdf"), width = 10, height = 8)



gtf_transcripts

one_gene<-"OVCH2"
find_transcript_ids<-function(one_gene){
  mini_gtf_df<-gtf_transcripts%>%filter(gene_name==one_gene)
  mini_gtf_df$transcript_id
}
list_of_transcripts<-find_transcript_ids(one_gene)
res_df_one_gene<-res_df%>%filter(transcript_id %in% list_of_transcripts )

# making bed file
gtf_transcripts_bed<-gtf_transcripts%>%filter(feature=="transcript") %>% select(seqname,start,end,transcript_id,score,strand)
dim(gtf_transcripts_bed)
colnames(gtf_transcripts_bed)

write.table(gtf_transcripts_bed,file = "gtf_transcripts_bed_output.bed",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
head(gtf_transcripts_bed)

colname1 <-c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'HG00268.merged.bam', 'HG00358.merged.bam', 'HG01352.merged.bam', 'HG01890.merged.bam', 'HG02059.merged.bam', 'HG02106.merged.bam', 'HG02282.merged.bam', 'HG02769.merged.bam', 'HG02818.merged.bam', 'HG03452.merged.bam', 'HG03456.merged.bam', 'HG03520.merged.bam', 'HG03807.merged.bam', 'HG04036.merged.bam', 'HG04217.merged.bam', 'NA19129.merged.bam', 'NA19434.merged.bam', 'NA19705.merged.bam', 'NA19836.merged.bam', 'NA20355.merged.bam', 'NA21487.merged.bam')

colname2 <- c("seqname", "start", "end", "name", "score", "strand")


nearbygenes <- read.csv("gtf_transcripts_bedwin_10kb.bed",sep="\t")

colnames(nearbygenes) <- c(colname1, colname2)

res_df_minexp <- res_df %>%
  dplyr::filter(!is.na(padj))%>%
  dplyr::filter(baseMean<100)


min_people<-read.csv("/scratch/Users/olde5615/data/graph21_RNA/RNA_hg38/sniffleswithnomorethen14peopleinonegenotype.csv")

# filter for sniffles near expressed genes and sniffles that vary between people with groups of no more than 14
nearbygenes_filter<-nearbygenes%>%filter(name %in% res_df_minexp$transcript_id)

#filters for genes and sniffle IDS
nearbygenes_filter_min_people<-nearbygenes_filter%>%filter(ID %in% min_people$row_0 )

n_distinct(nearbygenes_filter$name)
#[1] 73534
n_distinct(nearbygenes_filter$ID)
#[1] 33626

n_distinct(nearbygenes_filter_min_people$name)
n_distinct(nearbygenes_filter_min_people$ID)

write.csv(nearbygenes_filter_min_people,"candidate_sniffles.csv")





