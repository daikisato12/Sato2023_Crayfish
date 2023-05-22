library(tidyverse)
library(TCC)
setqwd("${data}/count/")

name="sen0_sen7"
#name="sen0_ton30"
#name="sen0_ton7"
#name="sen7_ton30"
#name="sen7_ton7"
#name="ton0_ton30"
#name="ton0_ton7"
#name="ton0_sen7"
#name="ton7_ton30"

input = paste0("gene_stringtie_non_novel_count_", name)
output = paste0("gene_stringtie_non_novel_", name, ".tsv")

group = c(1, 1, 1, 2, 2, 2)
iteration = 3
DE_method = "edger"
FDR = 0.05

count_data = read_tsv(input) %>% print()
count_matrix = count_data %>%
  tibble::column_to_rownames("gene_id") %>% 
  as.matrix()
tcc = new("TCC", count_matrix, group)
tcc = calcNormFactors(tcc, iteration = iteration)
tcc = estimateDE(tcc, test.method = DE_method, FDR = FDR)
result = getResult(tcc, sort = FALSE)
sorted = count_data %>%
  dplyr::left_join(result, by = "gene_id") %>%
  dplyr::arrange(rank) %>%
  print()
write_tsv(sorted, output)   #output file name
plot(tcc, FDR = FDR)
View (sorted)


