library(stringr)
library(purrr)
library(tidyr)
library(dplyr)
library(DESeq2)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
a_sniffle = args[1]
# sniffle with nearby.sig.gene
#a_sniffle="Sniffles2.INS.BAMA"
print(a_sniffle)

inDir <- "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/perSample/"
outdir<- "/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/DESeq2_results_a_sniffle/"

# Combine all samples into one count matrix
counts_df <- read.csv("/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/RNAseq_counts_allpeople.csv",row.names=1)
# remove NA19320 - not in graph genome
counts_df <- counts_df %>% select(-NA19320)

nearbygenes <- read.table("/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts/gtf_transcripts_bedwin_10kb.bed",sep="\t",header=FALSE, stringsAsFactors=FALSE)
# Define VCF header columns (first 9 columns)
vcf_cols <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')

# Define sample columns (excluding NA19320)
sample_cols <- c('HG00268.merged.bam', 'HG00358.merged.bam', 'HG01352.merged.bam', 
                 'HG01890.merged.bam', 'HG02059.merged.bam', 'HG02106.merged.bam', 
                 'HG02282.merged.bam', 'HG02769.merged.bam', 'HG02818.merged.bam', 
                 'HG03452.merged.bam', 'HG03456.merged.bam', 'HG03520.merged.bam', 
                 'HG03807.merged.bam', 'HG04036.merged.bam', 'HG04217.merged.bam', 
                 'NA19129.merged.bam', 'NA19434.merged.bam', 'NA19705.merged.bam', 
                 'NA19836.merged.bam', 'NA20355.merged.bam', 'NA21487.merged.bam')

# Define BED columns (last 6 columns)
bed_cols <- c("seqname", "start", "end", "name", "score", "strand")

# Combine all column names
all_cols <- c(vcf_cols, sample_cols, bed_cols)
print("test col n")
print(ncol(nearbygenes))
names(nearbygenes) <- all_cols

# Format metadata
meta <- read.csv("/scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/21_sample_metadata.csv")

# Filter for your favorite Sniffle (you can change this ID)
onesniffle <- nearbygenes %>% filter(ID==a_sniffle) %>% select(all_of(c(sample_cols))) %>% distinct()
print("onesniffle info")
print(colnames(onesniffle))
print(dim(onesniffle))
onesniffleT <- t(onesniffle)
head(onesniffle)
print("test col n 2")
print(ncol(onesniffleT))
colnames(onesniffleT) <- c("TE_genotype")
onesniffleT <- as.data.frame(onesniffleT)  
onesniffleT <- onesniffleT %>%
  mutate(TE_group = ifelse(TE_genotype == "0|0", 0L, 1L))

onesniffleT$bamFilename <- rownames(onesniffleT)
onesniffleT <- onesniffleT %>% separate_wider_delim(cols = bamFilename, delim = ".", names = c("SampleID", "merge", "bam")) 
onesniffleT <- onesniffleT %>% select(TE_group, SampleID)

meta2 <- merge(meta, onesniffleT, by = "SampleID")

print("TE_group_table")
print(table(meta2$TE_group))

# Filter metadata to only samples that exist in counts
meta3 <- as.data.frame(meta2) %>% filter(SampleID %in% colnames(counts_df))

# Prepare count matrix for DESeq2
rownames(counts_df) <- counts_df$Geneid
counts_df <- counts_df %>% select(-Geneid)
counts_df <- counts_df %>% select(all_of(meta3$SampleID)) 
head(counts_df)

# Set TE_group as factor
meta3$TE_group <- as.factor(meta3$TE_group)



# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_df, 
                              colData = meta3, 
                              design = ~ TE_group)

# Run DESeq2
DEdds <- DESeq(dds)

# Set your comparison groups (replace with your actual group names)
# Example: if TE_group has levels "0" and "1"
sample1 <- "0"  # Control/Reference group
sample2 <- "1"  # Test group

# Get results
res <- results(DEdds, contrast = c("TE_group", sample1, sample2))
# Filter your results
res_df <- as.data.frame(res)
res_df$transcript_id <- rownames(res_df)

transcripts_near_sniffle<-nearbygenes%>%filter(ID==a_sniffle)

# significant genes near sniffle - for individuals
res_sig<-res_df%>%filter(padj<0.1)
n_sig_transcripts_near_sniffle<-res_sig%>%filter(transcript_id %in% transcripts_near_sniffle$name)
if (nrow(n_sig_transcripts_near_sniffle)>0 ){
  outfile<-paste0(outdir,a_sniffle,".nearby.sig_genes.csv")
  write.csv(n_sig_transcripts_near_sniffle,outfile,row.names=FALSE) 
}

dim(res_sig)

res_sig_up<-res_sig%>%filter(log2FoldChange>0)
res_sig_down<-res_sig%>%filter(log2FoldChange<0)

changed_transcripts<-rbind(nrow(res_sig_up),nrow(res_sig_down))
rownames(changed_transcripts)<-c("up","down")
colnames(changed_transcripts)<-c(a_sniffle)

# save individual sniffle files
outfile<-paste0(outdir,a_sniffle,".n.sig_genes.csv")
write.csv(changed_transcripts,outfile,row.names = FALSE)



