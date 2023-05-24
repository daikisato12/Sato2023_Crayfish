#### load libraries ####
library(tidyverse)

#### Table 1 ####
##### get list of genes CAFE analysis ####
df_cafe_p <- read.delim("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/Gamma_branch_probabilities.tab", header = TRUE) %>%
  dplyr::rename(GeneID = X.FamilyID) %>%
  dplyr::select(!c(X.14., X)) %>%
  pivot_longer(cols = !GeneID, names_to = "node", values_to = "pval")
df_cafe_change <- read.delim("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/Gamma_change.tab", header = TRUE) %>%
  dplyr::rename(GeneID = FamilyID) %>%
  dplyr::select(!c(X.14.)) %>%
  pivot_longer(cols = !GeneID, names_to = "node", values_to = "change")
df_cafe <- inner_join(df_cafe_p %>%
                        filter(pval < 0.05), 
                      df_cafe_change) %>%
  mutate(sign = ifelse(change < 0, "-", ifelse(change > 0, "+", 0))) %>%
  dplyr::select(!change) %>%
  group_by(node, sign) %>%
  summarize(n = n())


df_cafe_Pcla_exp <- inner_join(df_cafe_p %>%
                                 filter(pval < 0.05), 
                               df_cafe_change) %>%
  mutate(sign = ifelse(change < 0, "-", ifelse(change > 0, "+", 0))) %>%
  filter(node == "Pcla.5.", sign == "+")


##### get list of genes DEG analysis ####
# list of Sendai-specific DEGs
list_degs_sendai <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval))

# list of Sapporo-specific DEGs
list_degs_sapporo <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         SapporoDay0_SapporoDay7_qval < 0.05 | SapporoDay0_SapporoDay31_qval  < 0.05)

list_degs_DEGs <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05 |
           SapporoDay0_SapporoDay7_qval < 0.05 | 
           SapporoDay0_SapporoDay31_qval  < 0.05 |
           SapporoDay7_SapporoDay31_qval  < 0.05)

df_cafe_num <- read.table("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/Gamma_count.tab", header = TRUE) %>%
  dplyr::rename(GeneID = FamilyID) %>%
  pivot_longer(cols = !GeneID, names_to = "node", values_to = "num")


list_degs_Pcla_exp <- bind_rows(list_degs_sendai, list_degs_sapporo)
# write.table(list_degs_sendai_exp, "../../ザリガニ論文_Sato_Makino/figures/tmp/Table1_sendai_exp.tsv", append = F, quote = F, row.names = F, sep = "\t")

list_degs_Pcla_exp2 <- list_degs_DEGs %>% #list_degs_Pcla_exp %>%
  dplyr::select(FlyID, GeneID, GeneName, AnnotatedName) %>%
  group_by(FlyID, GeneName) %>%
  summarize(DEG = n(),
            CrayfishID = paste(GeneID, collapse = ", ")) %>%
  inner_join(df_cafe_Pcla_exp %>%
               dplyr::rename(FlyID = GeneID)) %>%
  inner_join(df_cafe_num %>%
               dplyr::rename(FlyID = GeneID)) %>%
  dplyr::select(FlyID, GeneName, DEG, change, num, CrayfishID) %>%
  dplyr::rename(IncreaseNumGene = change,
                NumGene = num) %>%
  filter(DEG > 4)

write.table(list_degs_Pcla_exp2, 
            "../figures/tmp/Table1.tsv", 
            append = F, quote = F, row.names = F, sep = "\t")



