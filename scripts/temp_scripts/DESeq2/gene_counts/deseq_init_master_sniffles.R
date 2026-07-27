library(tidyverse)
# call in sniffles individually for the R array
mini_can_df<-read.csv("/scratch/Users/olde5615/data/pangenome21_phased/21_graphs_31MAY26/phased_vcf/candidate_sniffles_linear_genes.csv")


# creates full candidate output csv
output_df<-mini_can_df %>% select(ID) %>% distinct

write.csv(output_df,paste0("/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts_gene_id/", "just_candidate_sniffles_linear_genes.csv"),row.names=FALSE, quote=FALSE)


# creates mini candidate output csv
mini_outputdf <- head(output_df)

write.csv(mini_outputdf,paste0("/scratch/Users/olde5615/data/graph21_RNA/featureCounts_transcripts_gene_id/", "mini_just_candidate_sniffles_linear_genes.csv"),row.names=FALSE, quote=FALSE)
