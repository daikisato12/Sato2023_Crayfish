#### load libraries ####
library(tidyverse)

#### Table S2 ####
##### load dataset ####
# KEGG result of Sapporo-specific DEGs
df_kegg_sapporo <- read.delim("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/KEGG_enrichment_SapporoSpecificDEGs.tsv", header = TRUE, sep = "\t") %>%
  mutate(Type = "KEGG",
         Aim = "Sapporo-specific") %>%
  arrange(qvalue)
# KEGG result of Sendai-specific DEGs
df_kegg_sendai <- read.delim("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/KEGG_enrichment_SendaiSpecificDEGs.tsv", header = TRUE, sep = "\t") %>%
  mutate(Type = "KEGG",
         Aim = "Sendai-specific") %>%
  arrange(qvalue)
# KEGG result of Shared DEGs
df_kegg_shared <- read.delim("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/KEGG_enrichment_SharedDEGs.tsv", header = TRUE, sep = "\t") %>%
  mutate(Type = "KEGG",
         Aim = "Shared") %>%
  arrange(qvalue)


# Reactome result of Sapporo-specific DEGs
df_reactome_sapporo <- read.delim("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/Reactome_enrichment_SapporoSpecificDEGs.tsv", header = TRUE, sep = "\t") %>%
  rename(genename = geneID) %>%
  mutate(Type = "Reactome",
         Aim = "Sapporo-specific") %>%
  arrange(qvalue)
# Reactome result of Sendai-specific DEGs
df_reactome_sendai <- read.delim("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/Reactome_enrichment_SendaiSpecificDEGs.tsv", header = TRUE, sep = "\t") %>%
  rename(genename = geneID) %>%
  mutate(Type = "Reactome",
         Aim = "Sendai-specific") %>%
  arrange(qvalue)
# Reactome result of Shared DEGs
df_reactome_shared <- read.delim("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/Reactome_enrichment_SharedDEGs.tsv", header = TRUE, sep = "\t") %>%
  rename(genename = geneID) %>%
  mutate(Type = "Reactome",
         Aim = "Shared") %>%
  arrange(qvalue)

df_TableS2 <- bind_rows(df_kegg_sapporo, df_kegg_sendai) %>%
  bind_rows(df_kegg_shared) %>%
  bind_rows(df_reactome_sapporo) %>%
  bind_rows(df_reactome_sendai) %>%
  bind_rows(df_reactome_shared) %>%
  relocate(Type, Aim, everything()) %>%
  dplyr::select(!geneID)
  

write.table(df_TableS2, 
            "../figures/tmp/TableS2.tsv", 
            append = F, quote = F, row.names = F, sep = "\t")
