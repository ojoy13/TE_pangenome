library(stringr)
library(purrr)
library(tidyr)
library(dplyr)
library(DESeq2)
library(ggplot2)

inDir <- "/TE_pangenome/outputs/featureCounts/perSample/"

# Get list of ALL featureCounts files
files_list <- list.files(path = inDir, pattern = "\\.txt$", full.names = FALSE)

# Create list to store all sample data
df_list <- list()

# Loop through ALL files (no single-sample preprocessing needed)
for (i in seq_along(files_list)) {   
  fn <- paste0(inDir, files_list[i])
  fc <- read.table(fn, header = TRUE, skip = 1)
  
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
meta <- read.csv("/TE_pangenome/outputs/DESEq2/21_sample_metadata.csv")
sniffle <- read.csv("/TE_pangenome/outputs/DESEq2/deseq2_metadata.csv")

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
#counts_df <- counts_df %>% select(-Geneid)
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
res <- results(DEdds, contrast = c("TE_group", "0", "1"))
head(res)

# Convert to dataframe
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Load GTF for gene annotations
gtf <- "/path/to/dir/gencode.v38.ComprehensiveAnnotation.Main.gtf"

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

# Filter for gene features
gtf_genes <- gtf_df %>% filter(feature == "gene")

# Extract gene_id from the attribute column
gtf_genes <- gtf_genes %>%
  mutate(Geneid = str_extract(attribute, 'gene_id "([^"]+)"', group = 1))

# making bed file - CONVERT to 0-based coordinates
gtf_geneID_bed <- gtf_genes %>%
  filter(feature == "gene") %>% 
  mutate(start = start - 1) %>%  # Convert to 0-based for BED
  select(seqname, start, end, Geneid, score, strand)

# Check the first few rows
head(gtf_geneID_bed)
dim(gtf_geneID_bed)

write.table(gtf_geneID_bed, file = "gtf_geneID_bed_output.bed", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Run bedtools window from bash (you need to run this separately or with system())
# system("bedtools window -w 10000 -a /scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/pangenome_21_phased_nosymbolic.vcf.gz -b gtf_geneID_bed_output.bed > gtf_genes_bedwin_10kb.bed")

# Then read the results
nearbygenes <- read.csv("/TE_pangenome/outputs/DESEq2/gtf_genes_bedwin_10kb.bed", sep="\t")

colname1 <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 
              'HG00268.merged.bam', 'HG00358.merged.bam', 'HG01352.merged.bam', 
              'HG01890.merged.bam', 'HG02059.merged.bam', 'HG02106.merged.bam', 
              'HG02282.merged.bam', 'HG02769.merged.bam', 'HG02818.merged.bam', 
              'HG03452.merged.bam', 'HG03456.merged.bam', 'HG03520.merged.bam', 
              'HG03807.merged.bam', 'HG04036.merged.bam', 'HG04217.merged.bam', 
              'NA19129.merged.bam', 'NA19434.merged.bam', 'NA19705.merged.bam', 
              'NA19836.merged.bam', 'NA20355.merged.bam', 'NA21487.merged.bam') 

colname2 <- c("seqname", "start", "end", "name", "score", "strand")
colnames(nearbygenes) <- c(colname1, colname2)

# Check matching
cat("Number of genes in nearbygenes:", n_distinct(nearbygenes$name), "\n")
cat("Number of genes in res_df:", n_distinct(res_df$gene_id), "\n")
cat("Genes that match:", sum(nearbygenes$name %in% res_df$gene_id), "\n")

# Filter for HIGH expression genes
res_df_highexp <- res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::filter(baseMean > 100)

cat("Number of high expression genes:", nrow(res_df_highexp), "\n")

# Read the file with sniffles that have no more than 14 people in one genotype
min_people <- read.csv("/TE_pangenome/outputs/DESEq2/sniffleswithnomorethen14peopleinonegenotype.csv")

# Filter for sniffles near expressed genes
nearbygenes_filter <- nearbygenes %>% 
  filter(name %in% res_df_highexp$gene_id)

# Filters for genes and sniffle IDs
nearbygenes_filter_min_people <- nearbygenes_filter %>% 
  filter(ID %in% min_people$row_0)

# Count distinct genes
cat("Distinct genes in nearbygenes_filter:", n_distinct(nearbygenes_filter$name), "\n")
cat("Distinct sniffles in nearbygenes_filter:", n_distinct(nearbygenes_filter$ID), "\n")
cat("Distinct genes after min_people filter:", n_distinct(nearbygenes_filter_min_people$name), "\n")
cat("Distinct sniffles after min_people filter:", n_distinct(nearbygenes_filter_min_people$ID), "\n")

write.csv(nearbygenes_filter_min_people, "candidate_sniffles_linear_genes.csv", row.names = FALSE)
