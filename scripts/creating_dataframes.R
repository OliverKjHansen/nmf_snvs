library("callr")
library("tidyverse")


args<-commandArgs(trailingOnly = TRUE) 


background_counts <- args[1]
snv_counts  <- args[2]


path_to_output_snv <- args[3]
# path_to_input <- args[4]




background_matrix <- read.delim(background_counts, sep="\t", header = FALSE)
colnames(background_matrix) <- c("region","kmer","counts")
snv_matrix <- read.delim(snv_counts, sep="\t", header = FALSE)
colnames(snv_matrix) <- c("region","kmer","counts")

ms_bck_matrix <- pivot_wider(background_matrix, names_from = kmer, values_from = counts, values_fill = 0, names_sort = TRUE)
ms_bck_matrix <- ms_bck_matrix[str_order(ms_bck_matrix$region, numeric = TRUE),]

ms_snv_matrix <- pivot_wider(snv_matrix, names_from = kmer, values_from = counts, values_fill = 0, names_sort = TRUE)
ms_snv_matrix <- ms_snv_matrix[str_order(ms_snv_matrix$region, numeric = TRUE),]


ms_snv_matrix <- ms_snv_matrix %>% remove_rownames %>% column_to_rownames(var="region")
ms_bck_matrix <- ms_bck_matrix %>% remove_rownames %>% column_to_rownames(var="region")

ms_snv_df <- list("count_matrix" = ms_snv_matrix,
             "background_matrix" = ms_bck_matrix)
  


saveRDS(ms_snv_df, path_to_output_del)
