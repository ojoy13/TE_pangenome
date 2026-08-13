library(tidyverse)
# call in sniffles individually for the R array
mini_can_df<-read.csv("/TE_pangenome/outputs/DESEq2/candidate_sniffles.csv")


# creates full candidate output csv
output_df<-mini_can_df %>% select(ID) %>% distinct

write.csv(output_df,paste0("/TE_pangenome/outputs/DESEq2/", "just_candidate_sniffles.csv"),row.names=FALSE, quote=FALSE)


# creates mini candidate output csv
mini_outputdf <- head(output_df)

write.csv(mini_outputdf,paste0("/TE_pangenome/outputs/DESEq2/", "mini_just_candidate_sniffles.csv"),row.names=FALSE, quote=FALSE)
