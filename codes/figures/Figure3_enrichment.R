#### load libraries ####
# rm(list = ls(all.names = TRUE))
# library(tidyverse)
# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)

organism = "org.Dm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# BiocManager::install("enrichplot", force = TRUE)
# install.packages("ggupset")
# install.packages("ggnewscale")
library(ggupset)
library(ggnewscale)
library(enrichplot)

#### load genes ####
setwd("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/codes/")

# reading dataset
genelist_background = read.delim("../data/blastp/Allgenes_LG_crayfish_newID_flyID.txt", h=F)
genelist_background <- as.vector(genelist_background[,1]) %>% unique()

# reading dataset
genelist_rnaseq = read.delim("../data/RNA-seq/unique_mapping/rawdata/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104_flyID.txt", h=F)
genelist_rnaseq <- as.vector(genelist_rnaseq[,1]) %>% unique()

genelist_DEG2 = read.delim("../data/RNA-seq/unique_mapping/DEGs/RNAseq_DEG2_CrayfishID_flyID.txt", header=F)
genelist_DEG2 <- as.vector(genelist_DEG2[,1]) %>% unique()

genelist_PBS = read.delim("../data/PBS/PBS_mincov4_maxcov16_1kbp-top1_2combsharedloci_sharedgene_flyID.txt", header=F)
genelist_PBS <- as.vector(genelist_PBS[,1]) %>% unique()

genelist_PBS_DEG2_shared = read.delim("../data/PBS/1kbp/PBS_mincov4_maxcov16_1kbp-top1_2combsharedloci_DEGs_sharedgene_flyID.txt", header=F)
genelist_PBS_DEG2_shared <- as.vector(genelist_PBS_DEG2_shared[,1]) %>% unique()


#### PBS ####
##### BP #####
go_enrich_BP_PBS <- enrichGO(gene = genelist_PBS,
                             universe = genelist_background,
                             OrgDb = organism, 
                             keyType = 'FLYBASE',
                             readable = T,
                             ont = "BP",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.10)

df_go_enrich_BP_PBS <- go_enrich_BP_PBS@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_BP_PBS, "../data/enrichment/clusterProfiler_PBS/GO_BP_PBS.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_PBS <- upsetplot(go_enrich_BP_PBS)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/UpsetPlot_BP.pdf", g1_BP_PBS, h=6, w=10)

g2_BP_PBS <- barplot(go_enrich_BP_PBS, 
                     drop = TRUE, 
                     showCategory = 10, 
                     title = "GO Biological Pathways",
                     font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/BarPlot_BP.pdf", g2_BP_PBS, h=3, w=5)

g3_BP_PBS <- dotplot(go_enrich_BP_PBS, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/DotPlot_BP.pdf", g3_BP_PBS, h=4, w=6)

g4_BP_PBS <- emapplot(go_enrich_BP_PBS %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/EncrichmentMap_BP.pdf", g4_BP_PBS, h=6, w=8)

g5_BP_PBS <- goplot(go_enrich_BP_PBS, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/EncrichedGOGraph_BP.pdf", g5_BP_PBS, h=6, w=8)

g6_BP_PBS <- cnetplot(go_enrich_BP_PBS, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/CategoryNet_BP.pdf", g6_BP_PBS, h=12, w=10)

##### CC #####
go_enrich_CC_PBS <- enrichGO(gene = genelist_PBS,
                             universe = genelist_background,
                             OrgDb = organism, 
                             keyType = 'FLYBASE',
                             readable = T,
                             ont = "CC",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.10)

df_go_enrich_CC_PBS <- go_enrich_CC_PBS@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_CC_PBS, "../data/enrichment/clusterProfiler_PBS/GO_CC_PBS.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_PBS <- upsetplot(go_enrich_CC_PBS)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/UpsetPlot_CC.pdf", g1_CC_PBS, h=6, w=10)

g2_CC_PBS <- barplot(go_enrich_CC_PBS, 
                     drop = TRUE, 
                     showCategory = 10, 
                     title = "GO Cellular Components",
                     font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/BarPlot_CC.pdf", g2_CC_PBS, h=3, w=5)

g3_CC_PBS <- dotplot(go_enrich_CC_PBS, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/DotPlot_CC.pdf", g3_CC_PBS, h=4, w=5)

g4_CC_PBS <- emapplot(go_enrich_CC_PBS %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/EncrichmentMap_CC.pdf", g4_CC_PBS, h=6, w=8)

g5_CC_PBS <- goplot(go_enrich_CC_PBS, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/EncrichedGOGraph_CC.pdf", g5_CC_PBS, h=6, w=8)

g6_CC_PBS <- cnetplot(go_enrich_CC_PBS, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/CategoryNet_CC.pdf", g6_CC_PBS, h=12, w=10)

##### MF #####
go_enrich_MF_PBS <- enrichGO(gene = genelist_PBS,
                             universe = genelist_background,
                             OrgDb = organism, 
                             keyType = 'FLYBASE',
                             readable = T,
                             ont = "MF",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.10)

df_go_enrich_MF_PBS <- go_enrich_MF_PBS@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_MF_PBS, "../data/enrichment/clusterProfiler_PBS/GO_MF_PBS.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_PBS <- upsetplot(go_enrich_MF_PBS)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/UpsetPlot_MF.pdf", g1_MF_PBS, h=6, w=10)

g2_MF_PBS <- barplot(go_enrich_MF_PBS, 
                     drop = TRUE, 
                     showCategory = 10, 
                     title = "GO Molecular Functions",
                     font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/BarPlot_MF.pdf", g2_MF_PBS, h=3, w=5)

g3_MF_PBS <- dotplot(go_enrich_MF_PBS, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/DotPlot_MF.pdf", g3_MF_PBS, h=4, w=6)

g4_MF_PBS <- emapplot(go_enrich_MF_PBS %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/EncrichmentMap_MF.pdf", g4_MF_PBS, h=6, w=8)

g5_MF_PBS <- goplot(go_enrich_MF_PBS, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/EncrichedGOGraph_MF.pdf", g5_MF_PBS, h=6, w=8)

g6_MF_PBS <- cnetplot(go_enrich_MF_PBS, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS/CategoryNet_MF.pdf", g6_MF_PBS, h=12, w=10)


#### DEG2 ####
##### BP #####
go_enrich_BP_DEG2 <- enrichGO(gene = genelist_DEG2,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "BP",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_BP_DEG2 <- go_enrich_BP_DEG2@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_DEG2, "../data/enrichment/clusterProfiler_DEG2/GO_BP_DEG2.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_DEG2 <- upsetplot(go_enrich_BP_DEG2)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/UpsetPlot_BP.pdf", g1_BP_DEG2, h=6, w=10)

g2_BP_DEG2 <- barplot(go_enrich_BP_DEG2, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Biological Pathways",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/BarPlot_BP.pdf", g2_BP_DEG2, h=3, w=5)

g3_BP_DEG2 <- dotplot(go_enrich_BP_DEG2, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/DotPlot_BP.pdf", g3_BP_DEG2, h=4, w=6)

g4_BP_DEG2 <- emapplot(go_enrich_BP_DEG2 %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/EncrichmentMap_BP.pdf", g4_BP_DEG2, h=6, w=8)

g5_BP_DEG2 <- goplot(go_enrich_BP_DEG2, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/EncrichedGOGraph_BP.pdf", g5_BP_DEG2, h=6, w=8)

g6_BP_DEG2 <- cnetplot(go_enrich_BP_DEG2, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/CategoryNet_BP.pdf", g6_BP_DEG2, h=12, w=10)

##### CC #####
go_enrich_CC_DEG2 <- enrichGO(gene = genelist_DEG2,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "CC",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_CC_DEG2 <- go_enrich_CC_DEG2@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_DEG2, "../data/enrichment/clusterProfiler_DEG2/GO_CC_DEG2.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_DEG2 <- upsetplot(go_enrich_CC_DEG2)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/UpsetPlot_CC.pdf", g1_CC_DEG2, h=6, w=10)

g2_CC_DEG2 <- barplot(go_enrich_CC_DEG2, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Cellular Components",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/BarPlot_CC.pdf", g2_CC_DEG2, h=3, w=5)

g3_CC_DEG2 <- dotplot(go_enrich_CC_DEG2, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/DotPlot_CC.pdf", g3_CC_DEG2, h=4, w=5)

g4_CC_DEG2 <- emapplot(go_enrich_CC_DEG2 %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/EncrichmentMap_CC.pdf", g4_CC_DEG2, h=6, w=8)

g5_CC_DEG2 <- goplot(go_enrich_CC_DEG2, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/EncrichedGOGraph_CC.pdf", g5_CC_DEG2, h=6, w=8)

g6_CC_DEG2 <- cnetplot(go_enrich_CC_DEG2, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/CategoryNet_CC.pdf", g6_CC_DEG2, h=12, w=10)

##### MF #####
go_enrich_MF_DEG2 <- enrichGO(gene = genelist_DEG2,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "MF",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_MF_DEG2 <- go_enrich_MF_DEG2@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_DEG2, "../data/enrichment/clusterProfiler_DEG2/GO_MF_DEG2.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_DEG2 <- upsetplot(go_enrich_MF_DEG2)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/UpsetPlot_MF.pdf", g1_MF_DEG2, h=6, w=10)

g2_MF_DEG2 <- barplot(go_enrich_MF_DEG2, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Molecular Functions",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/BarPlot_MF.pdf", g2_MF_DEG2, h=3, w=5)

g3_MF_DEG2 <- dotplot(go_enrich_MF_DEG2, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/DotPlot_MF.pdf", g3_MF_DEG2, h=4, w=6)

g4_MF_DEG2 <- emapplot(go_enrich_MF_DEG2 %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/EncrichmentMap_MF.pdf", g4_MF_DEG2, h=6, w=8)

g5_MF_DEG2 <- goplot(go_enrich_MF_DEG2, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/EncrichedGOGraph_MF.pdf", g5_MF_DEG2, h=6, w=8)

g6_MF_DEG2 <- cnetplot(go_enrich_MF_DEG2, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_DEG2/CategoryNet_MF.pdf", g6_MF_DEG2, h=12, w=10)



#### PBS_DEG2_shared ####
##### BP #####
go_enrich_BP_PBS_DEG2_shared <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                         universe = genelist_background,
                                         OrgDb = organism, 
                                         keyType = 'FLYBASE',
                                         readable = T,
                                         ont = "BP",
                                         pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.10)

df_go_enrich_BP_PBS_DEG2_shared <- go_enrich_BP_PBS_DEG2_shared@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_BP_PBS_DEG2_shared, "../data/enrichment/clusterProfiler_PBS_DEG2_shared/GO_BP_PBS_DEG2_shared.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_PBS_DEG2_shared <- upsetplot(go_enrich_BP_PBS_DEG2_shared)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/UpsetPlot_BP.pdf", g1_BP_PBS_DEG2_shared, h=6, w=10)

g2_BP_PBS_DEG2_shared <- barplot(go_enrich_BP_PBS_DEG2_shared, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Biological Pathways",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/BarPlot_BP.pdf", g2_BP_PBS_DEG2_shared, h=3, w=5)

g3_BP_PBS_DEG2_shared <- dotplot(go_enrich_BP_PBS_DEG2_shared, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/DotPlot_BP.pdf", g3_BP_PBS_DEG2_shared, h=4, w=6)

g4_BP_PBS_DEG2_shared <- emapplot(go_enrich_BP_PBS_DEG2_shared %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/EncrichmentMap_BP.pdf", g4_BP_PBS_DEG2_shared, h=6, w=8)

g5_BP_PBS_DEG2_shared <- goplot(go_enrich_BP_PBS_DEG2_shared, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/EncrichedGOGraph_BP.pdf", g5_BP_PBS_DEG2_shared, h=6, w=8)

g6_BP_PBS_DEG2_shared <- cnetplot(go_enrich_BP_PBS_DEG2_shared, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/CategoryNet_BP.pdf", g6_BP_PBS_DEG2_shared, h=12, w=10)

##### CC #####
go_enrich_CC_PBS_DEG2_shared <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                         universe = genelist_background,
                                         OrgDb = organism, 
                                         keyType = 'FLYBASE',
                                         readable = T,
                                         ont = "CC",
                                         pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.10)

df_go_enrich_CC_PBS_DEG2_shared <- go_enrich_CC_PBS_DEG2_shared@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_CC_PBS_DEG2_shared, "../data/enrichment/clusterProfiler_PBS_DEG2_shared/GO_CC_PBS_DEG2_shared.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_PBS_DEG2_shared <- upsetplot(go_enrich_CC_PBS_DEG2_shared)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/UpsetPlot_CC.pdf", g1_CC_PBS_DEG2_shared, h=6, w=10)

g2_CC_PBS_DEG2_shared <- barplot(go_enrich_CC_PBS_DEG2_shared, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Cellular Components",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/BarPlot_CC.pdf", g2_CC_PBS_DEG2_shared, h=3, w=5)

g3_CC_PBS_DEG2_shared <- dotplot(go_enrich_CC_PBS_DEG2_shared, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/DotPlot_CC.pdf", g3_CC_PBS_DEG2_shared, h=4, w=5)

g4_CC_PBS_DEG2_shared <- emapplot(go_enrich_CC_PBS_DEG2_shared %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/EncrichmentMap_CC.pdf", g4_CC_PBS_DEG2_shared, h=6, w=8)

g5_CC_PBS_DEG2_shared <- goplot(go_enrich_CC_PBS_DEG2_shared, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/EncrichedGOGraph_CC.pdf", g5_CC_PBS_DEG2_shared, h=6, w=8)

g6_CC_PBS_DEG2_shared <- cnetplot(go_enrich_CC_PBS_DEG2_shared, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/CategoryNet_CC.pdf", g6_CC_PBS_DEG2_shared, h=12, w=10)

##### MF #####
go_enrich_MF_PBS_DEG2_shared <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                         universe = genelist_background,
                                         OrgDb = organism, 
                                         keyType = 'FLYBASE',
                                         readable = T,
                                         ont = "MF",
                                         pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.10)

df_go_enrich_MF_PBS_DEG2_shared <- go_enrich_MF_PBS_DEG2_shared@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_MF_PBS_DEG2_shared, "../data/enrichment/clusterProfiler_PBS_DEG2_shared/GO_MF_PBS_DEG2_shared.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_PBS_DEG2_shared <- upsetplot(go_enrich_MF_PBS_DEG2_shared)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/UpsetPlot_MF.pdf", g1_MF_PBS_DEG2_shared, h=6, w=10)

g2_MF_PBS_DEG2_shared <- barplot(go_enrich_MF_PBS_DEG2_shared, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Molecular Functions",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/BarPlot_MF.pdf", g2_MF_PBS_DEG2_shared, h=3, w=5)

g3_MF_PBS_DEG2_shared <- dotplot(go_enrich_MF_PBS_DEG2_shared, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/DotPlot_MF.pdf", g3_MF_PBS_DEG2_shared, h=4, w=6)

g4_MF_PBS_DEG2_shared <- emapplot(go_enrich_MF_PBS_DEG2_shared %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/EncrichmentMap_MF.pdf", g4_MF_PBS_DEG2_shared, h=6, w=8)

g5_MF_PBS_DEG2_shared <- goplot(go_enrich_MF_PBS_DEG2_shared, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/EncrichedGOGraph_MF.pdf", g5_MF_PBS_DEG2_shared, h=6, w=8)

g6_MF_PBS_DEG2_shared <- cnetplot(go_enrich_MF_PBS_DEG2_shared, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared/CategoryNet_MF.pdf", g6_MF_PBS_DEG2_shared, h=12, w=10)


#### PBS_DEG2_shared RNA-seq background####
##### BP #####
go_enrich_BP_PBS_DEG2_shared_RNAseqback <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "BP",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_BP_PBS_DEG2_shared_RNAseqback <- go_enrich_BP_PBS_DEG2_shared_RNAseqback@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_BP_PBS_DEG2_shared_RNAseqback, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/GO_BP_PBS_DEG2_shared_RNAseqback.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_PBS_DEG2_shared <- upsetplot(go_enrich_BP_PBS_DEG2_shared_RNAseqback)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/UpsetPlot_BP.pdf", g1_BP_PBS_DEG2_shared, h=6, w=10)

g2_BP_PBS_DEG2_shared <- barplot(go_enrich_BP_PBS_DEG2_shared_RNAseqback, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Biological Pathways",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/BarPlot_BP.pdf", g2_BP_PBS_DEG2_shared, h=3, w=5)

g3_BP_PBS_DEG2_shared <- dotplot(go_enrich_BP_PBS_DEG2_shared_RNAseqback, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/DotPlot_BP.pdf", g3_BP_PBS_DEG2_shared, h=4, w=6)

g4_BP_PBS_DEG2_shared <- emapplot(go_enrich_BP_PBS_DEG2_shared_RNAseqback %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/EncrichmentMap_BP.pdf", g4_BP_PBS_DEG2_shared, h=6, w=8)

g5_BP_PBS_DEG2_shared <- goplot(go_enrich_BP_PBS_DEG2_shared_RNAseqback, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/EncrichedGOGraph_BP.pdf", g5_BP_PBS_DEG2_shared, h=6, w=8)

g6_BP_PBS_DEG2_shared <- cnetplot(go_enrich_BP_PBS_DEG2_shared_RNAseqback, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/CategoryNet_BP.pdf", g6_BP_PBS_DEG2_shared, h=12, w=10)

##### CC #####
go_enrich_CC_PBS_DEG2_shared_RNAseqback <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "CC",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_CC_PBS_DEG2_shared_RNAseqback <- go_enrich_CC_PBS_DEG2_shared_RNAseqback@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_CC_PBS_DEG2_shared_RNAseqback, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/GO_CC_PBS_DEG2_shared_RNAseqback.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_PBS_DEG2_shared <- upsetplot(go_enrich_CC_PBS_DEG2_shared_RNAseqback)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/UpsetPlot_CC.pdf", g1_CC_PBS_DEG2_shared, h=6, w=10)

g2_CC_PBS_DEG2_shared <- barplot(go_enrich_CC_PBS_DEG2_shared_RNAseqback, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Cellular Components",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/BarPlot_CC.pdf", g2_CC_PBS_DEG2_shared, h=3, w=5)

g3_CC_PBS_DEG2_shared <- dotplot(go_enrich_CC_PBS_DEG2_shared_RNAseqback, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/DotPlot_CC.pdf", g3_CC_PBS_DEG2_shared, h=4, w=5)

g4_CC_PBS_DEG2_shared <- emapplot(go_enrich_CC_PBS_DEG2_shared_RNAseqback %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/EncrichmentMap_CC.pdf", g4_CC_PBS_DEG2_shared, h=6, w=8)

g5_CC_PBS_DEG2_shared <- goplot(go_enrich_CC_PBS_DEG2_shared_RNAseqback, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/EncrichedGOGraph_CC.pdf", g5_CC_PBS_DEG2_shared, h=6, w=8)

g6_CC_PBS_DEG2_shared <- cnetplot(go_enrich_CC_PBS_DEG2_shared_RNAseqback, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/CategoryNet_CC.pdf", g6_CC_PBS_DEG2_shared, h=12, w=10)

##### MF #####
go_enrich_MF_PBS_DEG2_shared_RNAseqback <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "MF",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_MF_PBS_DEG2_shared_RNAseqback <- go_enrich_MF_PBS_DEG2_shared_RNAseqback@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_MF_PBS_DEG2_shared_RNAseqback, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/GO_MF_PBS_DEG2_shared_RNAseqback.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_PBS_DEG2_shared <- upsetplot(go_enrich_MF_PBS_DEG2_shared_RNAseqback)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/UpsetPlot_MF.pdf", g1_MF_PBS_DEG2_shared, h=6, w=10)

g2_MF_PBS_DEG2_shared <- barplot(go_enrich_MF_PBS_DEG2_shared_RNAseqback, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Molecular Functions",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/BarPlot_MF.pdf", g2_MF_PBS_DEG2_shared, h=3, w=5)

g3_MF_PBS_DEG2_shared <- dotplot(go_enrich_MF_PBS_DEG2_shared_RNAseqback, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/DotPlot_MF.pdf", g3_MF_PBS_DEG2_shared, h=4, w=6)

g4_MF_PBS_DEG2_shared <- emapplot(go_enrich_MF_PBS_DEG2_shared_RNAseqback %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/EncrichmentMap_MF.pdf", g4_MF_PBS_DEG2_shared, h=6, w=8)

g5_MF_PBS_DEG2_shared <- goplot(go_enrich_MF_PBS_DEG2_shared_RNAseqback, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/EncrichedGOGraph_MF.pdf", g5_MF_PBS_DEG2_shared, h=6, w=8)

g6_MF_PBS_DEG2_shared <- cnetplot(go_enrich_MF_PBS_DEG2_shared_RNAseqback, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_RNAseq_background/CategoryNet_MF.pdf", g6_MF_PBS_DEG2_shared, h=12, w=10)


#### PBS_DEG2_shared DEG2 background####
##### BP #####
go_enrich_BP_PBS_DEG2_shared_DEG2back <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                  universe = genelist_DEG2,
                                                  OrgDb = organism, 
                                                  keyType = 'FLYBASE',
                                                  readable = T,
                                                  ont = "BP",
                                                  pvalueCutoff = 0.05, 
                                                  qvalueCutoff = 0.10)

df_go_enrich_BP_PBS_DEG2_shared_DEG2back <- go_enrich_BP_PBS_DEG2_shared_DEG2back@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_BP_PBS_DEG2_shared_DEG2back, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/GO_BP_PBS_DEG2_shared_DEG2back.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_PBS_DEG2_shared <- upsetplot(go_enrich_BP_PBS_DEG2_shared_DEG2back)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/UpsetPlot_BP.pdf", g1_BP_PBS_DEG2_shared, h=6, w=10)

g2_BP_PBS_DEG2_shared <- barplot(go_enrich_BP_PBS_DEG2_shared_DEG2back, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Biological Pathways",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/BarPlot_BP.pdf", g2_BP_PBS_DEG2_shared, h=3, w=5)

g3_BP_PBS_DEG2_shared <- dotplot(go_enrich_BP_PBS_DEG2_shared_DEG2back, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/DotPlot_BP.pdf", g3_BP_PBS_DEG2_shared, h=4, w=6)

g4_BP_PBS_DEG2_shared <- emapplot(go_enrich_BP_PBS_DEG2_shared_DEG2back %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/EncrichmentMap_BP.pdf", g4_BP_PBS_DEG2_shared, h=6, w=8)

g5_BP_PBS_DEG2_shared <- goplot(go_enrich_BP_PBS_DEG2_shared_DEG2back, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/EncrichedGOGraph_BP.pdf", g5_BP_PBS_DEG2_shared, h=6, w=8)

g6_BP_PBS_DEG2_shared <- cnetplot(go_enrich_BP_PBS_DEG2_shared_DEG2back, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/CategoryNet_BP.pdf", g6_BP_PBS_DEG2_shared, h=12, w=10)

##### CC #####
go_enrich_CC_PBS_DEG2_shared_DEG2back <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                  universe = genelist_DEG2,
                                                  OrgDb = organism, 
                                                  keyType = 'FLYBASE',
                                                  readable = T,
                                                  ont = "CC",
                                                  pvalueCutoff = 0.05, 
                                                  qvalueCutoff = 0.10)

df_go_enrich_CC_PBS_DEG2_shared_DEG2back <- go_enrich_CC_PBS_DEG2_shared_DEG2back@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_CC_PBS_DEG2_shared_DEG2back, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/GO_CC_PBS_DEG2_shared_DEG2back.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_PBS_DEG2_shared <- upsetplot(go_enrich_CC_PBS_DEG2_shared_DEG2back)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/UpsetPlot_CC.pdf", g1_CC_PBS_DEG2_shared, h=6, w=10)

g2_CC_PBS_DEG2_shared <- barplot(go_enrich_CC_PBS_DEG2_shared_DEG2back, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Cellular Components",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/BarPlot_CC.pdf", g2_CC_PBS_DEG2_shared, h=3, w=5)

g3_CC_PBS_DEG2_shared <- dotplot(go_enrich_CC_PBS_DEG2_shared_DEG2back, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/DotPlot_CC.pdf", g3_CC_PBS_DEG2_shared, h=4, w=5)

g4_CC_PBS_DEG2_shared <- emapplot(go_enrich_CC_PBS_DEG2_shared_DEG2back %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/EncrichmentMap_CC.pdf", g4_CC_PBS_DEG2_shared, h=6, w=8)

g5_CC_PBS_DEG2_shared <- goplot(go_enrich_CC_PBS_DEG2_shared_DEG2back, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/EncrichedGOGraph_CC.pdf", g5_CC_PBS_DEG2_shared, h=6, w=8)

g6_CC_PBS_DEG2_shared <- cnetplot(go_enrich_CC_PBS_DEG2_shared_DEG2back, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/CategoryNet_CC.pdf", g6_CC_PBS_DEG2_shared, h=12, w=10)

##### MF #####
go_enrich_MF_PBS_DEG2_shared_DEG2back <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                  universe = genelist_DEG2,
                                                  OrgDb = organism, 
                                                  keyType = 'FLYBASE',
                                                  readable = T,
                                                  ont = "MF",
                                                  pvalueCutoff = 0.05, 
                                                  qvalueCutoff = 0.10)

df_go_enrich_MF_PBS_DEG2_shared_DEG2back <- go_enrich_MF_PBS_DEG2_shared_DEG2back@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_MF_PBS_DEG2_shared_DEG2back, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/GO_MF_PBS_DEG2_shared_DEG2back.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_PBS_DEG2_shared <- upsetplot(go_enrich_MF_PBS_DEG2_shared_DEG2back)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/UpsetPlot_MF.pdf", g1_MF_PBS_DEG2_shared, h=6, w=10)

g2_MF_PBS_DEG2_shared <- barplot(go_enrich_MF_PBS_DEG2_shared_DEG2back, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Molecular Functions",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/BarPlot_MF.pdf", g2_MF_PBS_DEG2_shared, h=3, w=5)

g3_MF_PBS_DEG2_shared <- dotplot(go_enrich_MF_PBS_DEG2_shared_DEG2back, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/DotPlot_MF.pdf", g3_MF_PBS_DEG2_shared, h=4, w=6)

g4_MF_PBS_DEG2_shared <- emapplot(go_enrich_MF_PBS_DEG2_shared_DEG2back %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/EncrichmentMap_MF.pdf", g4_MF_PBS_DEG2_shared, h=6, w=8)

g5_MF_PBS_DEG2_shared <- goplot(go_enrich_MF_PBS_DEG2_shared_DEG2back, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/EncrichedGOGraph_MF.pdf", g5_MF_PBS_DEG2_shared, h=6, w=8)

g6_MF_PBS_DEG2_shared <- cnetplot(go_enrich_MF_PBS_DEG2_shared_DEG2back, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_DEG2_background/CategoryNet_MF.pdf", g6_MF_PBS_DEG2_shared, h=12, w=10)


#### PBS_DEG2_shared PBS background####
##### BP #####
go_enrich_BP_PBS_DEG2_shared_PBSback <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                 universe = genelist_PBS,
                                                 OrgDb = organism, 
                                                 keyType = 'FLYBASE',
                                                 readable = T,
                                                 ont = "BP",
                                                 pvalueCutoff = 0.05, 
                                                 qvalueCutoff = 0.10)

df_go_enrich_BP_PBS_DEG2_shared_PBSback <- go_enrich_BP_PBS_DEG2_shared_PBSback@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_BP_PBS_DEG2_shared_PBSback, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/GO_BP_PBS_DEG2_shared_PBSback.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_PBS_DEG2_shared <- upsetplot(go_enrich_BP_PBS_DEG2_shared_PBSback)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/UpsetPlot_BP.pdf", g1_BP_PBS_DEG2_shared, h=6, w=10)

g2_BP_PBS_DEG2_shared <- barplot(go_enrich_BP_PBS_DEG2_shared_PBSback, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Biological Pathways",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/BarPlot_BP.pdf", g2_BP_PBS_DEG2_shared, h=3, w=5)

g3_BP_PBS_DEG2_shared <- dotplot(go_enrich_BP_PBS_DEG2_shared_PBSback, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/DotPlot_BP.pdf", g3_BP_PBS_DEG2_shared, h=4, w=6)

g4_BP_PBS_DEG2_shared <- emapplot(go_enrich_BP_PBS_DEG2_shared_PBSback %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/EncrichmentMap_BP.pdf", g4_BP_PBS_DEG2_shared, h=6, w=8)

g5_BP_PBS_DEG2_shared <- goplot(go_enrich_BP_PBS_DEG2_shared_PBSback, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/EncrichedGOGraph_BP.pdf", g5_BP_PBS_DEG2_shared, h=6, w=8)

g6_BP_PBS_DEG2_shared <- cnetplot(go_enrich_BP_PBS_DEG2_shared_PBSback, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/CategoryNet_BP.pdf", g6_BP_PBS_DEG2_shared, h=12, w=10)

##### CC #####
go_enrich_CC_PBS_DEG2_shared_PBSback <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                 universe = genelist_PBS,
                                                 OrgDb = organism, 
                                                 keyType = 'FLYBASE',
                                                 readable = T,
                                                 ont = "CC",
                                                 pvalueCutoff = 0.05, 
                                                 qvalueCutoff = 0.10)

df_go_enrich_CC_PBS_DEG2_shared_PBSback <- go_enrich_CC_PBS_DEG2_shared_PBSback@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_CC_PBS_DEG2_shared_PBSback, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/GO_CC_PBS_DEG2_shared_PBSback.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_PBS_DEG2_shared <- upsetplot(go_enrich_CC_PBS_DEG2_shared_PBSback)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/UpsetPlot_CC.pdf", g1_CC_PBS_DEG2_shared, h=6, w=10)

g2_CC_PBS_DEG2_shared <- barplot(go_enrich_CC_PBS_DEG2_shared_PBSback, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Cellular Components",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/BarPlot_CC.pdf", g2_CC_PBS_DEG2_shared, h=3, w=5)

g3_CC_PBS_DEG2_shared <- dotplot(go_enrich_CC_PBS_DEG2_shared_PBSback, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/DotPlot_CC.pdf", g3_CC_PBS_DEG2_shared, h=4, w=5)

g4_CC_PBS_DEG2_shared <- emapplot(go_enrich_CC_PBS_DEG2_shared_PBSback %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/EncrichmentMap_CC.pdf", g4_CC_PBS_DEG2_shared, h=6, w=8)

g5_CC_PBS_DEG2_shared <- goplot(go_enrich_CC_PBS_DEG2_shared_PBSback, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/EncrichedGOGraph_CC.pdf", g5_CC_PBS_DEG2_shared, h=6, w=8)

g6_CC_PBS_DEG2_shared <- cnetplot(go_enrich_CC_PBS_DEG2_shared_PBSback, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/CategoryNet_CC.pdf", g6_CC_PBS_DEG2_shared, h=12, w=10)

##### MF #####
go_enrich_MF_PBS_DEG2_shared_PBSback <- enrichGO(gene = genelist_PBS_DEG2_shared,
                                                 universe = genelist_PBS,
                                                 OrgDb = organism, 
                                                 keyType = 'FLYBASE',
                                                 readable = T,
                                                 ont = "MF",
                                                 pvalueCutoff = 0.05, 
                                                 qvalueCutoff = 0.10)

df_go_enrich_MF_PBS_DEG2_shared_PBSback <- go_enrich_MF_PBS_DEG2_shared_PBSback@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_go_enrich_MF_PBS_DEG2_shared_PBSback, "../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/GO_MF_PBS_DEG2_shared_PBSback.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_PBS_DEG2_shared <- upsetplot(go_enrich_MF_PBS_DEG2_shared_PBSback)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/UpsetPlot_MF.pdf", g1_MF_PBS_DEG2_shared, h=6, w=10)

g2_MF_PBS_DEG2_shared <- barplot(go_enrich_MF_PBS_DEG2_shared_PBSback, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Molecular Functions",
                                 font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/BarPlot_MF.pdf", g2_MF_PBS_DEG2_shared, h=3, w=5)

g3_MF_PBS_DEG2_shared <- dotplot(go_enrich_MF_PBS_DEG2_shared_PBSback, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/DotPlot_MF.pdf", g3_MF_PBS_DEG2_shared, h=4, w=6)

g4_MF_PBS_DEG2_shared <- emapplot(go_enrich_MF_PBS_DEG2_shared_PBSback %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/EncrichmentMap_MF.pdf", g4_MF_PBS_DEG2_shared, h=6, w=8)

g5_MF_PBS_DEG2_shared <- goplot(go_enrich_MF_PBS_DEG2_shared_PBSback, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/EncrichedGOGraph_MF.pdf", g5_MF_PBS_DEG2_shared, h=6, w=8)

g6_MF_PBS_DEG2_shared <- cnetplot(go_enrich_MF_PBS_DEG2_shared_PBSback, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PBS_DEG2_shared_PBS_background/CategoryNet_MF.pdf", g6_MF_PBS_DEG2_shared, h=12, w=10)









#### Module ####
genelist_CTRG = read.delim("../figures/tmp/WGCNA_Module.tsv", header = T) 
genelist_CTRG <- as.vector(genelist_CTRG$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_CTRG <- enrichGO(gene = genelist_CTRG,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "BP",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_BP_CTRG <- go_enrich_BP_CTRG@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_CTRG, "../data/enrichment/clusterProfiler_CTRG/GO_BP_CTRG.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_CTRG <- upsetplot(go_enrich_BP_CTRG)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/UpsetPlot_BP.pdf", g1_BP_CTRG, h=6, w=10)

g2_BP_CTRG <- barplot(go_enrich_BP_CTRG, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Biological Pathways",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/BarPlot_BP.pdf", g2_BP_CTRG, h=3, w=5)

g3_BP_CTRG <- dotplot(go_enrich_BP_CTRG, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/DotPlot_BP.pdf", g3_BP_CTRG, h=4, w=6)

g4_BP_CTRG <- emapplot(go_enrich_BP_CTRG %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/EncrichmentMap_BP.pdf", g4_BP_CTRG, h=6, w=8)

g5_BP_CTRG <- goplot(go_enrich_BP_CTRG, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/EncrichedGOGraph_BP.pdf", g5_BP_CTRG, h=6, w=8)

g6_BP_CTRG <- cnetplot(go_enrich_BP_CTRG, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/CategoryNet_BP.pdf", g6_BP_CTRG, h=12, w=10)

##### CC #####
go_enrich_CC_CTRG <- enrichGO(gene = genelist_CTRG,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "CC",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_CC_CTRG <- go_enrich_CC_CTRG@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_CTRG, "../data/enrichment/clusterProfiler_CTRG/GO_CC_CTRG.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_CTRG <- upsetplot(go_enrich_CC_CTRG)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/UpsetPlot_CC.pdf", g1_CC_CTRG, h=6, w=10)

g2_CC_CTRG <- barplot(go_enrich_CC_CTRG, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Cellular Components",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/BarPlot_CC.pdf", g2_CC_CTRG, h=3, w=5)

g3_CC_CTRG <- dotplot(go_enrich_CC_CTRG, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/DotPlot_CC.pdf", g3_CC_CTRG, h=4, w=5)

g4_CC_CTRG <- emapplot(go_enrich_CC_CTRG %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/EncrichmentMap_CC.pdf", g4_CC_CTRG, h=6, w=8)

g5_CC_CTRG <- goplot(go_enrich_CC_CTRG, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/EncrichedGOGraph_CC.pdf", g5_CC_CTRG, h=6, w=8)

g6_CC_CTRG <- cnetplot(go_enrich_CC_CTRG, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/CategoryNet_CC.pdf", g6_CC_CTRG, h=12, w=10)

##### MF #####
go_enrich_MF_CTRG <- enrichGO(gene = genelist_CTRG,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "MF",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_MF_CTRG <- go_enrich_MF_CTRG@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_CTRG, "../data/enrichment/clusterProfiler_CTRG/GO_MF_CTRG.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_CTRG <- upsetplot(go_enrich_MF_CTRG)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/UpsetPlot_MF.pdf", g1_MF_CTRG, h=6, w=10)

g2_MF_CTRG <- barplot(go_enrich_MF_CTRG, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Molecular Functions",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/BarPlot_MF.pdf", g2_MF_CTRG, h=3, w=5)

g3_MF_CTRG <- dotplot(go_enrich_MF_CTRG, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/DotPlot_MF.pdf", g3_MF_CTRG, h=4, w=6)

g4_MF_CTRG <- emapplot(go_enrich_MF_CTRG %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/EncrichmentMap_MF.pdf", g4_MF_CTRG, h=6, w=8)

g5_MF_CTRG <- goplot(go_enrich_MF_CTRG, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/EncrichedGOGraph_MF.pdf", g5_MF_CTRG, h=6, w=8)

g6_MF_CTRG <- cnetplot(go_enrich_MF_CTRG, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_CTRG/CategoryNet_MF.pdf", g6_MF_CTRG, h=12, w=10)




#### Sendai specific DEGs ####
genelist_SendaiDEGs = read.delim("../figures/tmp/TableS5_SendaiSpecificDEGs.tsv", header = T) 
genelist_SendaiDEGs <- as.vector(genelist_SendaiDEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiDEGs <- enrichGO(gene = genelist_SendaiDEGs,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "BP",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiDEGs <- go_enrich_BP_SendaiDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiDEGs, "../data/enrichment/clusterProfiler_SendaiDEGs/GO_BP_SendaiDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiDEGs <- upsetplot(go_enrich_BP_SendaiDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/UpsetPlot_BP.pdf", g1_BP_SendaiDEGs, h=6, w=10)

g2_BP_SendaiDEGs <- barplot(go_enrich_BP_SendaiDEGs, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Biological Pathways",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/BarPlot_BP.pdf", g2_BP_SendaiDEGs, h=3, w=5)

g3_BP_SendaiDEGs <- dotplot(go_enrich_BP_SendaiDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/DotPlot_BP.pdf", g3_BP_SendaiDEGs, h=4, w=6)

g4_BP_SendaiDEGs <- emapplot(go_enrich_BP_SendaiDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiDEGs, h=6, w=8)

g5_BP_SendaiDEGs <- goplot(go_enrich_BP_SendaiDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiDEGs, h=6, w=8)

g6_BP_SendaiDEGs <- cnetplot(go_enrich_BP_SendaiDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/CategoryNet_BP.pdf", g6_BP_SendaiDEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiDEGs <- enrichGO(gene = genelist_SendaiDEGs,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "CC",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiDEGs <- go_enrich_CC_SendaiDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiDEGs, "../data/enrichment/clusterProfiler_SendaiDEGs/GO_CC_SendaiDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiDEGs <- upsetplot(go_enrich_CC_SendaiDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/UpsetPlot_CC.pdf", g1_CC_SendaiDEGs, h=6, w=10)

g2_CC_SendaiDEGs <- barplot(go_enrich_CC_SendaiDEGs, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Cellular Components",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/BarPlot_CC.pdf", g2_CC_SendaiDEGs, h=3, w=5)

g3_CC_SendaiDEGs <- dotplot(go_enrich_CC_SendaiDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/DotPlot_CC.pdf", g3_CC_SendaiDEGs, h=4, w=5)

g4_CC_SendaiDEGs <- emapplot(go_enrich_CC_SendaiDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiDEGs, h=6, w=8)

g5_CC_SendaiDEGs <- goplot(go_enrich_CC_SendaiDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiDEGs, h=6, w=8)

g6_CC_SendaiDEGs <- cnetplot(go_enrich_CC_SendaiDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/CategoryNet_CC.pdf", g6_CC_SendaiDEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiDEGs <- enrichGO(gene = genelist_SendaiDEGs,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "MF",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiDEGs <- go_enrich_MF_SendaiDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiDEGs, "../data/enrichment/clusterProfiler_SendaiDEGs/GO_MF_SendaiDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiDEGs <- upsetplot(go_enrich_MF_SendaiDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/UpsetPlot_MF.pdf", g1_MF_SendaiDEGs, h=6, w=10)

g2_MF_SendaiDEGs <- barplot(go_enrich_MF_SendaiDEGs, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Molecular Functions",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/BarPlot_MF.pdf", g2_MF_SendaiDEGs, h=3, w=5)

g3_MF_SendaiDEGs <- dotplot(go_enrich_MF_SendaiDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/DotPlot_MF.pdf", g3_MF_SendaiDEGs, h=4, w=6)

g4_MF_SendaiDEGs <- emapplot(go_enrich_MF_SendaiDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiDEGs, h=6, w=8)

g5_MF_SendaiDEGs <- goplot(go_enrich_MF_SendaiDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiDEGs, h=6, w=8)

g6_MF_SendaiDEGs <- cnetplot(go_enrich_MF_SendaiDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGs/CategoryNet_MF.pdf", g6_MF_SendaiDEGs, h=12, w=10)



#### Sapporo specific DEGs ####
genelist_SapporoDEGs = read.delim("../figures/tmp/TableS5_SapporoSpecificDEGs.tsv", header = T) 
genelist_SapporoDEGs <- as.vector(genelist_SapporoDEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SapporoDEGs <- enrichGO(gene = genelist_SapporoDEGs,
                                    universe = genelist_rnaseq,
                                    OrgDb = organism, 
                                    keyType = 'FLYBASE',
                                    readable = T,
                                    ont = "BP",
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SapporoDEGs <- go_enrich_BP_SapporoDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SapporoDEGs, "../data/enrichment/clusterProfiler_SapporoDEGs/GO_BP_SapporoDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SapporoDEGs <- upsetplot(go_enrich_BP_SapporoDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/UpsetPlot_BP.pdf", g1_BP_SapporoDEGs, h=6, w=10)

g2_BP_SapporoDEGs <- barplot(go_enrich_BP_SapporoDEGs, 
                            drop = TRUE, 
                            showCategory = 10, 
                            title = "GO Biological Pathways",
                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/BarPlot_BP.pdf", g2_BP_SapporoDEGs, h=3, w=5)

g3_BP_SapporoDEGs <- dotplot(go_enrich_BP_SapporoDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/DotPlot_BP.pdf", g3_BP_SapporoDEGs, h=4, w=6)

g4_BP_SapporoDEGs <- emapplot(go_enrich_BP_SapporoDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/EncrichmentMap_BP.pdf", g4_BP_SapporoDEGs, h=6, w=8)

g5_BP_SapporoDEGs <- goplot(go_enrich_BP_SapporoDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SapporoDEGs, h=6, w=8)

g6_BP_SapporoDEGs <- cnetplot(go_enrich_BP_SapporoDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/CategoryNet_BP.pdf", g6_BP_SapporoDEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SapporoDEGs <- enrichGO(gene = genelist_SapporoDEGs,
                                    universe = genelist_rnaseq,
                                    OrgDb = organism, 
                                    keyType = 'FLYBASE',
                                    readable = T,
                                    ont = "CC",
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SapporoDEGs <- go_enrich_CC_SapporoDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SapporoDEGs, "../data/enrichment/clusterProfiler_SapporoDEGs/GO_CC_SapporoDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SapporoDEGs <- upsetplot(go_enrich_CC_SapporoDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/UpsetPlot_CC.pdf", g1_CC_SapporoDEGs, h=6, w=10)

g2_CC_SapporoDEGs <- barplot(go_enrich_CC_SapporoDEGs, 
                            drop = TRUE, 
                            showCategory = 10, 
                            title = "GO Cellular Components",
                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/BarPlot_CC.pdf", g2_CC_SapporoDEGs, h=3, w=5)

g3_CC_SapporoDEGs <- dotplot(go_enrich_CC_SapporoDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/DotPlot_CC.pdf", g3_CC_SapporoDEGs, h=4, w=5)

g4_CC_SapporoDEGs <- emapplot(go_enrich_CC_SapporoDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/EncrichmentMap_CC.pdf", g4_CC_SapporoDEGs, h=6, w=8)

g5_CC_SapporoDEGs <- goplot(go_enrich_CC_SapporoDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SapporoDEGs, h=6, w=8)

g6_CC_SapporoDEGs <- cnetplot(go_enrich_CC_SapporoDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/CategoryNet_CC.pdf", g6_CC_SapporoDEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SapporoDEGs <- enrichGO(gene = genelist_SapporoDEGs,
                                    universe = genelist_rnaseq,
                                    OrgDb = organism, 
                                    keyType = 'FLYBASE',
                                    readable = T,
                                    ont = "MF",
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SapporoDEGs <- go_enrich_MF_SapporoDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SapporoDEGs, "../data/enrichment/clusterProfiler_SapporoDEGs/GO_MF_SapporoDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SapporoDEGs <- upsetplot(go_enrich_MF_SapporoDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/UpsetPlot_MF.pdf", g1_MF_SapporoDEGs, h=6, w=10)

g2_MF_SapporoDEGs <- barplot(go_enrich_MF_SapporoDEGs, 
                            drop = TRUE, 
                            showCategory = 10, 
                            title = "GO Molecular Functions",
                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/BarPlot_MF.pdf", g2_MF_SapporoDEGs, h=3, w=5)

g3_MF_SapporoDEGs <- dotplot(go_enrich_MF_SapporoDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/DotPlot_MF.pdf", g3_MF_SapporoDEGs, h=4, w=6)

g4_MF_SapporoDEGs <- emapplot(go_enrich_MF_SapporoDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/EncrichmentMap_MF.pdf", g4_MF_SapporoDEGs, h=6, w=8)

g5_MF_SapporoDEGs <- goplot(go_enrich_MF_SapporoDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SapporoDEGs, h=6, w=8)

g6_MF_SapporoDEGs <- cnetplot(go_enrich_MF_SapporoDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDEGs/CategoryNet_MF.pdf", g6_MF_SapporoDEGs, h=12, w=10)



#### SapporoDay31SpecificDEGs ####
genelist_SapporoDay31DEGs = read.delim("../figures/tmp/TableS5_SapporoDay31SpecificDEGs.tsv", header = T) 
genelist_SapporoDay31DEGs <- as.vector(genelist_SapporoDay31DEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SapporoDay31DEGs <- enrichGO(gene = genelist_SapporoDay31DEGs,
                                     universe = genelist_rnaseq,
                                     OrgDb = organism, 
                                     keyType = 'FLYBASE',
                                     readable = T,
                                     ont = "BP",
                                     pvalueCutoff = 0.05, 
                                     qvalueCutoff = 0.10)

df_go_enrich_BP_SapporoDay31DEGs <- go_enrich_BP_SapporoDay31DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SapporoDay31DEGs, "../data/enrichment/clusterProfiler_SapporoDay31DEGs/GO_BP_SapporoDay31DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SapporoDay31DEGs <- upsetplot(go_enrich_BP_SapporoDay31DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/UpsetPlot_BP.pdf", g1_BP_SapporoDay31DEGs, h=6, w=10)

g2_BP_SapporoDay31DEGs <- barplot(go_enrich_BP_SapporoDay31DEGs, 
                             drop = TRUE, 
                             showCategory = 10, 
                             title = "GO Biological Pathways",
                             font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/BarPlot_BP.pdf", g2_BP_SapporoDay31DEGs, h=3, w=5)

g3_BP_SapporoDay31DEGs <- dotplot(go_enrich_BP_SapporoDay31DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/DotPlot_BP.pdf", g3_BP_SapporoDay31DEGs, h=4, w=6)

g4_BP_SapporoDay31DEGs <- emapplot(go_enrich_BP_SapporoDay31DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/EncrichmentMap_BP.pdf", g4_BP_SapporoDay31DEGs, h=6, w=8)

g5_BP_SapporoDay31DEGs <- goplot(go_enrich_BP_SapporoDay31DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/EncrichedGOGraph_BP.pdf", g5_BP_SapporoDay31DEGs, h=6, w=8)

g6_BP_SapporoDay31DEGs <- cnetplot(go_enrich_BP_SapporoDay31DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/CategoryNet_BP.pdf", g6_BP_SapporoDay31DEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SapporoDay31DEGs <- enrichGO(gene = genelist_SapporoDay31DEGs,
                                     universe = genelist_rnaseq,
                                     OrgDb = organism, 
                                     keyType = 'FLYBASE',
                                     readable = T,
                                     ont = "CC",
                                     pvalueCutoff = 0.05, 
                                     qvalueCutoff = 0.10)

df_go_enrich_CC_SapporoDay31DEGs <- go_enrich_CC_SapporoDay31DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SapporoDay31DEGs, "../data/enrichment/clusterProfiler_SapporoDay31DEGs/GO_CC_SapporoDay31DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SapporoDay31DEGs <- upsetplot(go_enrich_CC_SapporoDay31DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/UpsetPlot_CC.pdf", g1_CC_SapporoDay31DEGs, h=6, w=10)

g2_CC_SapporoDay31DEGs <- barplot(go_enrich_CC_SapporoDay31DEGs, 
                             drop = TRUE, 
                             showCategory = 10, 
                             title = "GO Cellular Components",
                             font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/BarPlot_CC.pdf", g2_CC_SapporoDay31DEGs, h=3, w=5)

g3_CC_SapporoDay31DEGs <- dotplot(go_enrich_CC_SapporoDay31DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/DotPlot_CC.pdf", g3_CC_SapporoDay31DEGs, h=4, w=5)

g4_CC_SapporoDay31DEGs <- emapplot(go_enrich_CC_SapporoDay31DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/EncrichmentMap_CC.pdf", g4_CC_SapporoDay31DEGs, h=6, w=8)

g5_CC_SapporoDay31DEGs <- goplot(go_enrich_CC_SapporoDay31DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/EncrichedGOGraph_CC.pdf", g5_CC_SapporoDay31DEGs, h=6, w=8)

g6_CC_SapporoDay31DEGs <- cnetplot(go_enrich_CC_SapporoDay31DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/CategoryNet_CC.pdf", g6_CC_SapporoDay31DEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SapporoDay31DEGs <- enrichGO(gene = genelist_SapporoDay31DEGs,
                                     universe = genelist_rnaseq,
                                     OrgDb = organism, 
                                     keyType = 'FLYBASE',
                                     readable = T,
                                     ont = "MF",
                                     pvalueCutoff = 0.05, 
                                     qvalueCutoff = 0.10)

df_go_enrich_MF_SapporoDay31DEGs <- go_enrich_MF_SapporoDay31DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SapporoDay31DEGs, "../data/enrichment/clusterProfiler_SapporoDay31DEGs/GO_MF_SapporoDay31DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SapporoDay31DEGs <- upsetplot(go_enrich_MF_SapporoDay31DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/UpsetPlot_MF.pdf", g1_MF_SapporoDay31DEGs, h=6, w=10)

g2_MF_SapporoDay31DEGs <- barplot(go_enrich_MF_SapporoDay31DEGs, 
                             drop = TRUE, 
                             showCategory = 10, 
                             title = "GO Molecular Functions",
                             font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/BarPlot_MF.pdf", g2_MF_SapporoDay31DEGs, h=3, w=5)

g3_MF_SapporoDay31DEGs <- dotplot(go_enrich_MF_SapporoDay31DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/DotPlot_MF.pdf", g3_MF_SapporoDay31DEGs, h=4, w=6)

g4_MF_SapporoDay31DEGs <- emapplot(go_enrich_MF_SapporoDay31DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/EncrichmentMap_MF.pdf", g4_MF_SapporoDay31DEGs, h=6, w=8)

g5_MF_SapporoDay31DEGs <- goplot(go_enrich_MF_SapporoDay31DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/EncrichedGOGraph_MF.pdf", g5_MF_SapporoDay31DEGs, h=6, w=8)

g6_MF_SapporoDay31DEGs <- cnetplot(go_enrich_MF_SapporoDay31DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SapporoDay31DEGs/CategoryNet_MF.pdf", g6_MF_SapporoDay31DEGs, h=12, w=10)

#### SendaiOnlyDEGs ####
genelist_SendaiOnlyDEGs = read.delim("../figures/tmp/TableS5_SendaiOnlyDEGs.tsv", header = T) 
genelist_SendaiOnlyDEGs <- as.vector(genelist_SendaiOnlyDEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiOnlyDEGs <- enrichGO(gene = genelist_SendaiOnlyDEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "BP",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiOnlyDEGs <- go_enrich_BP_SendaiOnlyDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiOnlyDEGs, "../data/enrichment/clusterProfiler_SendaiOnlyDEGs/GO_BP_SendaiOnlyDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiOnlyDEGs <- upsetplot(go_enrich_BP_SendaiOnlyDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/UpsetPlot_BP.pdf", g1_BP_SendaiOnlyDEGs, h=6, w=10)

g2_BP_SendaiOnlyDEGs <- barplot(go_enrich_BP_SendaiOnlyDEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Biological Pathways",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/BarPlot_BP.pdf", g2_BP_SendaiOnlyDEGs, h=3, w=5)

g3_BP_SendaiOnlyDEGs <- dotplot(go_enrich_BP_SendaiOnlyDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/DotPlot_BP.pdf", g3_BP_SendaiOnlyDEGs, h=4, w=6)

g4_BP_SendaiOnlyDEGs <- emapplot(go_enrich_BP_SendaiOnlyDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiOnlyDEGs, h=6, w=8)

g5_BP_SendaiOnlyDEGs <- goplot(go_enrich_BP_SendaiOnlyDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiOnlyDEGs, h=6, w=8)

g6_BP_SendaiOnlyDEGs <- cnetplot(go_enrich_BP_SendaiOnlyDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/CategoryNet_BP.pdf", g6_BP_SendaiOnlyDEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiOnlyDEGs <- enrichGO(gene = genelist_SendaiOnlyDEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "CC",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiOnlyDEGs <- go_enrich_CC_SendaiOnlyDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiOnlyDEGs, "../data/enrichment/clusterProfiler_SendaiOnlyDEGs/GO_CC_SendaiOnlyDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiOnlyDEGs <- upsetplot(go_enrich_CC_SendaiOnlyDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/UpsetPlot_CC.pdf", g1_CC_SendaiOnlyDEGs, h=6, w=10)

g2_CC_SendaiOnlyDEGs <- barplot(go_enrich_CC_SendaiOnlyDEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Cellular Components",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/BarPlot_CC.pdf", g2_CC_SendaiOnlyDEGs, h=3, w=5)

g3_CC_SendaiOnlyDEGs <- dotplot(go_enrich_CC_SendaiOnlyDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/DotPlot_CC.pdf", g3_CC_SendaiOnlyDEGs, h=4, w=5)

g4_CC_SendaiOnlyDEGs <- emapplot(go_enrich_CC_SendaiOnlyDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiOnlyDEGs, h=6, w=8)

g5_CC_SendaiOnlyDEGs <- goplot(go_enrich_CC_SendaiOnlyDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiOnlyDEGs, h=6, w=8)

g6_CC_SendaiOnlyDEGs <- cnetplot(go_enrich_CC_SendaiOnlyDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/CategoryNet_CC.pdf", g6_CC_SendaiOnlyDEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiOnlyDEGs <- enrichGO(gene = genelist_SendaiOnlyDEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "MF",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiOnlyDEGs <- go_enrich_MF_SendaiOnlyDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiOnlyDEGs, "../data/enrichment/clusterProfiler_SendaiOnlyDEGs/GO_MF_SendaiOnlyDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiOnlyDEGs <- upsetplot(go_enrich_MF_SendaiOnlyDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/UpsetPlot_MF.pdf", g1_MF_SendaiOnlyDEGs, h=6, w=10)

g2_MF_SendaiOnlyDEGs <- barplot(go_enrich_MF_SendaiOnlyDEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Molecular Functions",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/BarPlot_MF.pdf", g2_MF_SendaiOnlyDEGs, h=3, w=5)

g3_MF_SendaiOnlyDEGs <- dotplot(go_enrich_MF_SendaiOnlyDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/DotPlot_MF.pdf", g3_MF_SendaiOnlyDEGs, h=4, w=6)

g4_MF_SendaiOnlyDEGs <- emapplot(go_enrich_MF_SendaiOnlyDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiOnlyDEGs, h=6, w=8)

g5_MF_SendaiOnlyDEGs <- goplot(go_enrich_MF_SendaiOnlyDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiOnlyDEGs, h=6, w=8)

g6_MF_SendaiOnlyDEGs <- cnetplot(go_enrich_MF_SendaiOnlyDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiOnlyDEGs/CategoryNet_MF.pdf", g6_MF_SendaiOnlyDEGs, h=12, w=10)



#### Population specific DEGs ####
genelist_PopulationSpecificDEGs = read.delim("../figures/tmp/TableS5_PopulationSpecificDEGs.tsv", header = T) 
genelist_PopulationSpecificDEGs <- as.vector(genelist_PopulationSpecificDEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_PopulationSpecificDEGs <- enrichGO(gene = genelist_PopulationSpecificDEGs,
                                     universe = genelist_rnaseq,
                                     OrgDb = organism, 
                                     keyType = 'FLYBASE',
                                     readable = T,
                                     ont = "BP",
                                     pvalueCutoff = 0.05, 
                                     qvalueCutoff = 0.10)

df_go_enrich_BP_PopulationSpecificDEGs <- go_enrich_BP_PopulationSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_PopulationSpecificDEGs, "../data/enrichment/clusterProfiler_PopulationSpecificDEGs/GO_BP_PopulationSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_PopulationSpecificDEGs <- upsetplot(go_enrich_BP_PopulationSpecificDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/UpsetPlot_BP.pdf", g1_BP_PopulationSpecificDEGs, h=6, w=10)

g2_BP_PopulationSpecificDEGs <- barplot(go_enrich_BP_PopulationSpecificDEGs, 
                             drop = TRUE, 
                             showCategory = 10, 
                             title = "GO Biological Pathways",
                             font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/BarPlot_BP.pdf", g2_BP_PopulationSpecificDEGs, h=3, w=5)

g3_BP_PopulationSpecificDEGs <- dotplot(go_enrich_BP_PopulationSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/DotPlot_BP.pdf", g3_BP_PopulationSpecificDEGs, h=4, w=6)

g4_BP_PopulationSpecificDEGs <- emapplot(go_enrich_BP_PopulationSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/EncrichmentMap_BP.pdf", g4_BP_PopulationSpecificDEGs, h=6, w=8)

g5_BP_PopulationSpecificDEGs <- goplot(go_enrich_BP_PopulationSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/EncrichedGOGraph_BP.pdf", g5_BP_PopulationSpecificDEGs, h=6, w=8)

g6_BP_PopulationSpecificDEGs <- cnetplot(go_enrich_BP_PopulationSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/CategoryNet_BP.pdf", g6_BP_PopulationSpecificDEGs, h=12, w=10)

##### CC #####
go_enrich_CC_PopulationSpecificDEGs <- enrichGO(gene = genelist_PopulationSpecificDEGs,
                                     universe = genelist_rnaseq,
                                     OrgDb = organism, 
                                     keyType = 'FLYBASE',
                                     readable = T,
                                     ont = "CC",
                                     pvalueCutoff = 0.05, 
                                     qvalueCutoff = 0.10)

df_go_enrich_CC_PopulationSpecificDEGs <- go_enrich_CC_PopulationSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_PopulationSpecificDEGs, "../data/enrichment/clusterProfiler_PopulationSpecificDEGs/GO_CC_PopulationSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_PopulationSpecificDEGs <- upsetplot(go_enrich_CC_PopulationSpecificDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/UpsetPlot_CC.pdf", g1_CC_PopulationSpecificDEGs, h=6, w=10)

g2_CC_PopulationSpecificDEGs <- barplot(go_enrich_CC_PopulationSpecificDEGs, 
                             drop = TRUE, 
                             showCategory = 10, 
                             title = "GO Cellular Components",
                             font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/BarPlot_CC.pdf", g2_CC_PopulationSpecificDEGs, h=3, w=5)

g3_CC_PopulationSpecificDEGs <- dotplot(go_enrich_CC_PopulationSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/DotPlot_CC.pdf", g3_CC_PopulationSpecificDEGs, h=4, w=5)

g4_CC_PopulationSpecificDEGs <- emapplot(go_enrich_CC_PopulationSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/EncrichmentMap_CC.pdf", g4_CC_PopulationSpecificDEGs, h=6, w=8)

g5_CC_PopulationSpecificDEGs <- goplot(go_enrich_CC_PopulationSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/EncrichedGOGraph_CC.pdf", g5_CC_PopulationSpecificDEGs, h=6, w=8)

g6_CC_PopulationSpecificDEGs <- cnetplot(go_enrich_CC_PopulationSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/CategoryNet_CC.pdf", g6_CC_PopulationSpecificDEGs, h=12, w=10)

##### MF #####
go_enrich_MF_PopulationSpecificDEGs <- enrichGO(gene = genelist_PopulationSpecificDEGs,
                                     universe = genelist_rnaseq,
                                     OrgDb = organism, 
                                     keyType = 'FLYBASE',
                                     readable = T,
                                     ont = "MF",
                                     pvalueCutoff = 0.05, 
                                     qvalueCutoff = 0.10)

df_go_enrich_MF_PopulationSpecificDEGs <- go_enrich_MF_PopulationSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_PopulationSpecificDEGs, "../data/enrichment/clusterProfiler_PopulationSpecificDEGs/GO_MF_PopulationSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_PopulationSpecificDEGs <- upsetplot(go_enrich_MF_PopulationSpecificDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/UpsetPlot_MF.pdf", g1_MF_PopulationSpecificDEGs, h=6, w=10)

g2_MF_PopulationSpecificDEGs <- barplot(go_enrich_MF_PopulationSpecificDEGs, 
                             drop = TRUE, 
                             showCategory = 10, 
                             title = "GO Molecular Functions",
                             font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/BarPlot_MF.pdf", g2_MF_PopulationSpecificDEGs, h=3, w=5)

g3_MF_PopulationSpecificDEGs <- dotplot(go_enrich_MF_PopulationSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/DotPlot_MF.pdf", g3_MF_PopulationSpecificDEGs, h=4, w=6)

g4_MF_PopulationSpecificDEGs <- emapplot(go_enrich_MF_PopulationSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/EncrichmentMap_MF.pdf", g4_MF_PopulationSpecificDEGs, h=6, w=8)

g5_MF_PopulationSpecificDEGs <- goplot(go_enrich_MF_PopulationSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/EncrichedGOGraph_MF.pdf", g5_MF_PopulationSpecificDEGs, h=6, w=8)

g6_MF_PopulationSpecificDEGs <- cnetplot(go_enrich_MF_PopulationSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopulationSpecificDEGs/CategoryNet_MF.pdf", g6_MF_PopulationSpecificDEGs, h=12, w=10)


#### SharedDEGs ####
genelist_SharedDEGs = read.delim("../figures/tmp/TableS5_SharedDEGs.tsv", header = T) 
genelist_SharedDEGs <- as.vector(genelist_SharedDEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SharedDEGs <- enrichGO(gene = genelist_SharedDEGs,
                                                                  universe = genelist_rnaseq,
                                                                  OrgDb = organism, 
                                                                  keyType = 'FLYBASE',
                                                                  readable = T,
                                                                  ont = "BP",
                                                                  pvalueCutoff = 0.05, 
                                                                  qvalueCutoff = 0.10)

df_go_enrich_BP_SharedDEGs <- go_enrich_BP_SharedDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SharedDEGs, "../data/enrichment/clusterProfiler_SharedDEGs/GO_BP_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SharedDEGs <- upsetplot(go_enrich_BP_SharedDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/UpsetPlot_BP.pdf", g1_BP_SharedDEGs, h=6, w=10)

g2_BP_SharedDEGs <- barplot(go_enrich_BP_SharedDEGs, 
                                                          drop = TRUE, 
                                                          showCategory = 10, 
                                                          title = "GO Biological Pathways",
                                                          font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/BarPlot_BP.pdf", g2_BP_SharedDEGs, h=3, w=5)

g3_BP_SharedDEGs <- dotplot(go_enrich_BP_SharedDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/DotPlot_BP.pdf", g3_BP_SharedDEGs, h=4, w=6)

g4_BP_SharedDEGs <- emapplot(go_enrich_BP_SharedDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/EncrichmentMap_BP.pdf", g4_BP_SharedDEGs, h=6, w=8)

g5_BP_SharedDEGs <- goplot(go_enrich_BP_SharedDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SharedDEGs, h=6, w=8)

g6_BP_SharedDEGs <- cnetplot(go_enrich_BP_SharedDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/CategoryNet_BP.pdf", g6_BP_SharedDEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SharedDEGs <- enrichGO(gene = genelist_SharedDEGs,
                                                                  universe = genelist_rnaseq,
                                                                  OrgDb = organism, 
                                                                  keyType = 'FLYBASE',
                                                                  readable = T,
                                                                  ont = "CC",
                                                                  pvalueCutoff = 0.05, 
                                                                  qvalueCutoff = 0.10)

df_go_enrich_CC_SharedDEGs <- go_enrich_CC_SharedDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SharedDEGs, "../data/enrichment/clusterProfiler_SharedDEGs/GO_CC_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SharedDEGs <- upsetplot(go_enrich_CC_SharedDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/UpsetPlot_CC.pdf", g1_CC_SharedDEGs, h=6, w=10)

g2_CC_SharedDEGs <- barplot(go_enrich_CC_SharedDEGs, 
                                                          drop = TRUE, 
                                                          showCategory = 10, 
                                                          title = "GO Cellular Components",
                                                          font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/BarPlot_CC.pdf", g2_CC_SharedDEGs, h=3, w=5)

g3_CC_SharedDEGs <- dotplot(go_enrich_CC_SharedDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/DotPlot_CC.pdf", g3_CC_SharedDEGs, h=4, w=5)

g4_CC_SharedDEGs <- emapplot(go_enrich_CC_SharedDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/EncrichmentMap_CC.pdf", g4_CC_SharedDEGs, h=6, w=8)

g5_CC_SharedDEGs <- goplot(go_enrich_CC_SharedDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SharedDEGs, h=6, w=8)

g6_CC_SharedDEGs <- cnetplot(go_enrich_CC_SharedDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/CategoryNet_CC.pdf", g6_CC_SharedDEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SharedDEGs <- enrichGO(gene = genelist_SharedDEGs,
                                                                  universe = genelist_rnaseq,
                                                                  OrgDb = organism, 
                                                                  keyType = 'FLYBASE',
                                                                  readable = T,
                                                                  ont = "MF",
                                                                  pvalueCutoff = 0.05, 
                                                                  qvalueCutoff = 0.10)

df_go_enrich_MF_SharedDEGs <- go_enrich_MF_SharedDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SharedDEGs, "../data/enrichment/clusterProfiler_SharedDEGs/GO_MF_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SharedDEGs <- upsetplot(go_enrich_MF_SharedDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/UpsetPlot_MF.pdf", g1_MF_SharedDEGs, h=6, w=10)

g2_MF_SharedDEGs <- barplot(go_enrich_MF_SharedDEGs, 
                                                          drop = TRUE, 
                                                          showCategory = 10, 
                                                          title = "GO Molecular Functions",
                                                          font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/BarPlot_MF.pdf", g2_MF_SharedDEGs, h=3, w=5)

g3_MF_SharedDEGs <- dotplot(go_enrich_MF_SharedDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/DotPlot_MF.pdf", g3_MF_SharedDEGs, h=4, w=6)

g4_MF_SharedDEGs <- emapplot(go_enrich_MF_SharedDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/EncrichmentMap_MF.pdf", g4_MF_SharedDEGs, h=6, w=8)

g5_MF_SharedDEGs <- goplot(go_enrich_MF_SharedDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SharedDEGs, h=6, w=8)

g6_MF_SharedDEGs <- cnetplot(go_enrich_MF_SharedDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SharedDEGs/CategoryNet_MF.pdf", g6_MF_SharedDEGs, h=12, w=10)




#### SendaiDEGsSapporo07NonDEGs ####
genelist_SendaiDEGsSapporo07NonDEGs = read.delim("../figures/tmp/TableS5_SendaiDEGsSapporo07NonDEGs.tsv", header = T) 
genelist_SendaiDEGsSapporo07NonDEGs <- as.vector(genelist_SendaiDEGsSapporo07NonDEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiDEGsSapporo07NonDEGs <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGs,
                                    universe = genelist_rnaseq,
                                    OrgDb = organism, 
                                    keyType = 'FLYBASE',
                                    readable = T,
                                    ont = "BP",
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiDEGsSapporo07NonDEGs <- go_enrich_BP_SendaiDEGsSapporo07NonDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiDEGsSapporo07NonDEGs, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/GO_BP_SendaiDEGsSapporo07NonDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiDEGsSapporo07NonDEGs <- upsetplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/UpsetPlot_BP.pdf", g1_BP_SendaiDEGsSapporo07NonDEGs, h=6, w=10)

g2_BP_SendaiDEGsSapporo07NonDEGs <- barplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs, 
                            drop = TRUE, 
                            showCategory = 10, 
                            title = "GO Biological Pathways",
                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/BarPlot_BP.pdf", g2_BP_SendaiDEGsSapporo07NonDEGs, h=3, w=5)

g3_BP_SendaiDEGsSapporo07NonDEGs <- dotplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/DotPlot_BP.pdf", g3_BP_SendaiDEGsSapporo07NonDEGs, h=4, w=6)

g4_BP_SendaiDEGsSapporo07NonDEGs <- emapplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiDEGsSapporo07NonDEGs, h=6, w=8)

g5_BP_SendaiDEGsSapporo07NonDEGs <- goplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiDEGsSapporo07NonDEGs, h=6, w=8)

g6_BP_SendaiDEGsSapporo07NonDEGs <- cnetplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/CategoryNet_BP.pdf", g6_BP_SendaiDEGsSapporo07NonDEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiDEGsSapporo07NonDEGs <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGs,
                                    universe = genelist_rnaseq,
                                    OrgDb = organism, 
                                    keyType = 'FLYBASE',
                                    readable = T,
                                    ont = "CC",
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiDEGsSapporo07NonDEGs <- go_enrich_CC_SendaiDEGsSapporo07NonDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiDEGsSapporo07NonDEGs, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/GO_CC_SendaiDEGsSapporo07NonDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiDEGsSapporo07NonDEGs <- upsetplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/UpsetPlot_CC.pdf", g1_CC_SendaiDEGsSapporo07NonDEGs, h=6, w=10)

g2_CC_SendaiDEGsSapporo07NonDEGs <- barplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs, 
                            drop = TRUE, 
                            showCategory = 10, 
                            title = "GO Cellular Components",
                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/BarPlot_CC.pdf", g2_CC_SendaiDEGsSapporo07NonDEGs, h=3, w=5)

g3_CC_SendaiDEGsSapporo07NonDEGs <- dotplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/DotPlot_CC.pdf", g3_CC_SendaiDEGsSapporo07NonDEGs, h=4, w=5)

g4_CC_SendaiDEGsSapporo07NonDEGs <- emapplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiDEGsSapporo07NonDEGs, h=6, w=8)

g5_CC_SendaiDEGsSapporo07NonDEGs <- goplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiDEGsSapporo07NonDEGs, h=6, w=8)

g6_CC_SendaiDEGsSapporo07NonDEGs <- cnetplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/CategoryNet_CC.pdf", g6_CC_SendaiDEGsSapporo07NonDEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiDEGsSapporo07NonDEGs <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGs,
                                    universe = genelist_rnaseq,
                                    OrgDb = organism, 
                                    keyType = 'FLYBASE',
                                    readable = T,
                                    ont = "MF",
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiDEGsSapporo07NonDEGs <- go_enrich_MF_SendaiDEGsSapporo07NonDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiDEGsSapporo07NonDEGs, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/GO_MF_SendaiDEGsSapporo07NonDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiDEGsSapporo07NonDEGs <- upsetplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/UpsetPlot_MF.pdf", g1_MF_SendaiDEGsSapporo07NonDEGs, h=6, w=10)

g2_MF_SendaiDEGsSapporo07NonDEGs <- barplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs, 
                            drop = TRUE, 
                            showCategory = 10, 
                            title = "GO Molecular Functions",
                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/BarPlot_MF.pdf", g2_MF_SendaiDEGsSapporo07NonDEGs, h=3, w=5)

g3_MF_SendaiDEGsSapporo07NonDEGs <- dotplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/DotPlot_MF.pdf", g3_MF_SendaiDEGsSapporo07NonDEGs, h=4, w=6)

g4_MF_SendaiDEGsSapporo07NonDEGs <- emapplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiDEGsSapporo07NonDEGs, h=6, w=8)

g5_MF_SendaiDEGsSapporo07NonDEGs <- goplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiDEGsSapporo07NonDEGs, h=6, w=8)

g6_MF_SendaiDEGsSapporo07NonDEGs <- cnetplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs/CategoryNet_MF.pdf", g6_MF_SendaiDEGsSapporo07NonDEGs, h=12, w=10)


#### SendaiDEGsorSapporo07NonDEGs031DEGs ####
genelist_SendaiDEGsorSapporo07NonDEGs031DEGs = read.delim("../figures/tmp/TableS5_SendaiDEGsorSapporo07NonDEGs031DEGs.tsv", header = T) 
genelist_SendaiDEGsorSapporo07NonDEGs031DEGs <- as.vector(genelist_SendaiDEGsorSapporo07NonDEGs031DEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- enrichGO(gene = genelist_SendaiDEGsorSapporo07NonDEGs031DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "BP",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, "../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/GO_BP_SendaiDEGsorSapporo07NonDEGs031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- upsetplot(go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/UpsetPlot_BP.pdf", g1_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=10)

g2_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- barplot(go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Biological Pathways",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/BarPlot_BP.pdf", g2_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, h=3, w=5)

g3_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- dotplot(go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/DotPlot_BP.pdf", g3_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, h=4, w=6)

g4_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- emapplot(go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=8)

g5_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- goplot(go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=8)

g6_BP_SendaiDEGsorSapporo07NonDEGs031DEGs <- cnetplot(go_enrich_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/CategoryNet_BP.pdf", g6_BP_SendaiDEGsorSapporo07NonDEGs031DEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- enrichGO(gene = genelist_SendaiDEGsorSapporo07NonDEGs031DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "CC",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, "../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/GO_CC_SendaiDEGsorSapporo07NonDEGs031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- upsetplot(go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/UpsetPlot_CC.pdf", g1_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=10)

g2_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- barplot(go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Cellular Components",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/BarPlot_CC.pdf", g2_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, h=3, w=5)

g3_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- dotplot(go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/DotPlot_CC.pdf", g3_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, h=4, w=5)

g4_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- emapplot(go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=8)

g5_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- goplot(go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=8)

g6_CC_SendaiDEGsorSapporo07NonDEGs031DEGs <- cnetplot(go_enrich_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/CategoryNet_CC.pdf", g6_CC_SendaiDEGsorSapporo07NonDEGs031DEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- enrichGO(gene = genelist_SendaiDEGsorSapporo07NonDEGs031DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "MF",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, "../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/GO_MF_SendaiDEGsorSapporo07NonDEGs031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- upsetplot(go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/UpsetPlot_MF.pdf", g1_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=10)

g2_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- barplot(go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Molecular Functions",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/BarPlot_MF.pdf", g2_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, h=3, w=5)

g3_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- dotplot(go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/DotPlot_MF.pdf", g3_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, h=4, w=6)

g4_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- emapplot(go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=8)

g5_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- goplot(go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, h=6, w=8)

g6_MF_SendaiDEGsorSapporo07NonDEGs031DEGs <- cnetplot(go_enrich_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsorSapporo07NonDEGs031DEGs/CategoryNet_MF.pdf", g6_MF_SendaiDEGsorSapporo07NonDEGs031DEGs, h=12, w=10)



#### SendaiNonDEGsSapporo07DEGs ####
genelist_SendaiNonDEGsSapporo07DEGs = read.delim("../figures/tmp/TableS5_SendaiNonDEGsSapporo07DEGs.tsv", header = T) 
genelist_SendaiNonDEGsSapporo07DEGs <- as.vector(genelist_SendaiNonDEGsSapporo07DEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiNonDEGsSapporo07DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo07DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "BP",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiNonDEGsSapporo07DEGs <- go_enrich_BP_SendaiNonDEGsSapporo07DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiNonDEGsSapporo07DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/GO_BP_SendaiNonDEGsSapporo07DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiNonDEGsSapporo07DEGs <- upsetplot(go_enrich_BP_SendaiNonDEGsSapporo07DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/UpsetPlot_BP.pdf", g1_BP_SendaiNonDEGsSapporo07DEGs, h=6, w=10)

g2_BP_SendaiNonDEGsSapporo07DEGs <- barplot(go_enrich_BP_SendaiNonDEGsSapporo07DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Biological Pathways",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/BarPlot_BP.pdf", g2_BP_SendaiNonDEGsSapporo07DEGs, h=3, w=5)

g3_BP_SendaiNonDEGsSapporo07DEGs <- dotplot(go_enrich_BP_SendaiNonDEGsSapporo07DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/DotPlot_BP.pdf", g3_BP_SendaiNonDEGsSapporo07DEGs, h=4, w=6)

g4_BP_SendaiNonDEGsSapporo07DEGs <- emapplot(go_enrich_BP_SendaiNonDEGsSapporo07DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiNonDEGsSapporo07DEGs, h=6, w=8)

g5_BP_SendaiNonDEGsSapporo07DEGs <- goplot(go_enrich_BP_SendaiNonDEGsSapporo07DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiNonDEGsSapporo07DEGs, h=6, w=8)

g6_BP_SendaiNonDEGsSapporo07DEGs <- cnetplot(go_enrich_BP_SendaiNonDEGsSapporo07DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/CategoryNet_BP.pdf", g6_BP_SendaiNonDEGsSapporo07DEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiNonDEGsSapporo07DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo07DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "CC",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiNonDEGsSapporo07DEGs <- go_enrich_CC_SendaiNonDEGsSapporo07DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiNonDEGsSapporo07DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/GO_CC_SendaiNonDEGsSapporo07DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiNonDEGsSapporo07DEGs <- upsetplot(go_enrich_CC_SendaiNonDEGsSapporo07DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/UpsetPlot_CC.pdf", g1_CC_SendaiNonDEGsSapporo07DEGs, h=6, w=10)

g2_CC_SendaiNonDEGsSapporo07DEGs <- barplot(go_enrich_CC_SendaiNonDEGsSapporo07DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Cellular Components",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/BarPlot_CC.pdf", g2_CC_SendaiNonDEGsSapporo07DEGs, h=3, w=5)

g3_CC_SendaiNonDEGsSapporo07DEGs <- dotplot(go_enrich_CC_SendaiNonDEGsSapporo07DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/DotPlot_CC.pdf", g3_CC_SendaiNonDEGsSapporo07DEGs, h=4, w=5)

g4_CC_SendaiNonDEGsSapporo07DEGs <- emapplot(go_enrich_CC_SendaiNonDEGsSapporo07DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiNonDEGsSapporo07DEGs, h=6, w=8)

g5_CC_SendaiNonDEGsSapporo07DEGs <- goplot(go_enrich_CC_SendaiNonDEGsSapporo07DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiNonDEGsSapporo07DEGs, h=6, w=8)

g6_CC_SendaiNonDEGsSapporo07DEGs <- cnetplot(go_enrich_CC_SendaiNonDEGsSapporo07DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/CategoryNet_CC.pdf", g6_CC_SendaiNonDEGsSapporo07DEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiNonDEGsSapporo07DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo07DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "MF",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiNonDEGsSapporo07DEGs <- go_enrich_MF_SendaiNonDEGsSapporo07DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiNonDEGsSapporo07DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/GO_MF_SendaiNonDEGsSapporo07DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiNonDEGsSapporo07DEGs <- upsetplot(go_enrich_MF_SendaiNonDEGsSapporo07DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/UpsetPlot_MF.pdf", g1_MF_SendaiNonDEGsSapporo07DEGs, h=6, w=10)

g2_MF_SendaiNonDEGsSapporo07DEGs <- barplot(go_enrich_MF_SendaiNonDEGsSapporo07DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Molecular Functions",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/BarPlot_MF.pdf", g2_MF_SendaiNonDEGsSapporo07DEGs, h=3, w=5)

g3_MF_SendaiNonDEGsSapporo07DEGs <- dotplot(go_enrich_MF_SendaiNonDEGsSapporo07DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/DotPlot_MF.pdf", g3_MF_SendaiNonDEGsSapporo07DEGs, h=4, w=6)

g4_MF_SendaiNonDEGsSapporo07DEGs <- emapplot(go_enrich_MF_SendaiNonDEGsSapporo07DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiNonDEGsSapporo07DEGs, h=6, w=8)

g5_MF_SendaiNonDEGsSapporo07DEGs <- goplot(go_enrich_MF_SendaiNonDEGsSapporo07DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiNonDEGsSapporo07DEGs, h=6, w=8)

g6_MF_SendaiNonDEGsSapporo07DEGs <- cnetplot(go_enrich_MF_SendaiNonDEGsSapporo07DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07DEGs/CategoryNet_MF.pdf", g6_MF_SendaiNonDEGsSapporo07DEGs, h=12, w=10)


#### SendaiNonDEGsSapporo07031DEGs ####
genelist_SendaiNonDEGsSapporo07031DEGs = read.delim("../figures/tmp/TableS5_SendaiNonDEGsSapporo07031DEGs.tsv", header = T) 
genelist_SendaiNonDEGsSapporo07031DEGs <- as.vector(genelist_SendaiNonDEGsSapporo07031DEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiNonDEGsSapporo07031DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo07031DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "BP",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiNonDEGsSapporo07031DEGs <- go_enrich_BP_SendaiNonDEGsSapporo07031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiNonDEGsSapporo07031DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/GO_BP_SendaiNonDEGsSapporo07031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiNonDEGsSapporo07031DEGs <- upsetplot(go_enrich_BP_SendaiNonDEGsSapporo07031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/UpsetPlot_BP.pdf", g1_BP_SendaiNonDEGsSapporo07031DEGs, h=6, w=10)

g2_BP_SendaiNonDEGsSapporo07031DEGs <- barplot(go_enrich_BP_SendaiNonDEGsSapporo07031DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Biological Pathways",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/BarPlot_BP.pdf", g2_BP_SendaiNonDEGsSapporo07031DEGs, h=3, w=5)

g3_BP_SendaiNonDEGsSapporo07031DEGs <- dotplot(go_enrich_BP_SendaiNonDEGsSapporo07031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/DotPlot_BP.pdf", g3_BP_SendaiNonDEGsSapporo07031DEGs, h=4, w=6)

g4_BP_SendaiNonDEGsSapporo07031DEGs <- emapplot(go_enrich_BP_SendaiNonDEGsSapporo07031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiNonDEGsSapporo07031DEGs, h=6, w=8)

g5_BP_SendaiNonDEGsSapporo07031DEGs <- goplot(go_enrich_BP_SendaiNonDEGsSapporo07031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiNonDEGsSapporo07031DEGs, h=6, w=8)

g6_BP_SendaiNonDEGsSapporo07031DEGs <- cnetplot(go_enrich_BP_SendaiNonDEGsSapporo07031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/CategoryNet_BP.pdf", g6_BP_SendaiNonDEGsSapporo07031DEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiNonDEGsSapporo07031DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo07031DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "CC",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiNonDEGsSapporo07031DEGs <- go_enrich_CC_SendaiNonDEGsSapporo07031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiNonDEGsSapporo07031DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/GO_CC_SendaiNonDEGsSapporo07031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiNonDEGsSapporo07031DEGs <- upsetplot(go_enrich_CC_SendaiNonDEGsSapporo07031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/UpsetPlot_CC.pdf", g1_CC_SendaiNonDEGsSapporo07031DEGs, h=6, w=10)

g2_CC_SendaiNonDEGsSapporo07031DEGs <- barplot(go_enrich_CC_SendaiNonDEGsSapporo07031DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Cellular Components",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/BarPlot_CC.pdf", g2_CC_SendaiNonDEGsSapporo07031DEGs, h=3, w=5)

g3_CC_SendaiNonDEGsSapporo07031DEGs <- dotplot(go_enrich_CC_SendaiNonDEGsSapporo07031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/DotPlot_CC.pdf", g3_CC_SendaiNonDEGsSapporo07031DEGs, h=4, w=5)

g4_CC_SendaiNonDEGsSapporo07031DEGs <- emapplot(go_enrich_CC_SendaiNonDEGsSapporo07031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiNonDEGsSapporo07031DEGs, h=6, w=8)

g5_CC_SendaiNonDEGsSapporo07031DEGs <- goplot(go_enrich_CC_SendaiNonDEGsSapporo07031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiNonDEGsSapporo07031DEGs, h=6, w=8)

g6_CC_SendaiNonDEGsSapporo07031DEGs <- cnetplot(go_enrich_CC_SendaiNonDEGsSapporo07031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/CategoryNet_CC.pdf", g6_CC_SendaiNonDEGsSapporo07031DEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiNonDEGsSapporo07031DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo07031DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "MF",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiNonDEGsSapporo07031DEGs <- go_enrich_MF_SendaiNonDEGsSapporo07031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiNonDEGsSapporo07031DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/GO_MF_SendaiNonDEGsSapporo07031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiNonDEGsSapporo07031DEGs <- upsetplot(go_enrich_MF_SendaiNonDEGsSapporo07031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/UpsetPlot_MF.pdf", g1_MF_SendaiNonDEGsSapporo07031DEGs, h=6, w=10)

g2_MF_SendaiNonDEGsSapporo07031DEGs <- barplot(go_enrich_MF_SendaiNonDEGsSapporo07031DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Molecular Functions",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/BarPlot_MF.pdf", g2_MF_SendaiNonDEGsSapporo07031DEGs, h=3, w=5)

g3_MF_SendaiNonDEGsSapporo07031DEGs <- dotplot(go_enrich_MF_SendaiNonDEGsSapporo07031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/DotPlot_MF.pdf", g3_MF_SendaiNonDEGsSapporo07031DEGs, h=4, w=6)

g4_MF_SendaiNonDEGsSapporo07031DEGs <- emapplot(go_enrich_MF_SendaiNonDEGsSapporo07031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiNonDEGsSapporo07031DEGs, h=6, w=8)

g5_MF_SendaiNonDEGsSapporo07031DEGs <- goplot(go_enrich_MF_SendaiNonDEGsSapporo07031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiNonDEGsSapporo07031DEGs, h=6, w=8)

g6_MF_SendaiNonDEGsSapporo07031DEGs <- cnetplot(go_enrich_MF_SendaiNonDEGsSapporo07031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo07031DEGs/CategoryNet_MF.pdf", g6_MF_SendaiNonDEGsSapporo07031DEGs, h=12, w=10)



#### SendaiNonDEGsSapporo031DEGs ####
genelist_SendaiNonDEGsSapporo031DEGs = read.delim("../figures/tmp/TableS5_SendaiNonDEGsSapporo031DEGs.tsv", header = T) 
genelist_SendaiNonDEGsSapporo031DEGs <- as.vector(genelist_SendaiNonDEGsSapporo031DEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiNonDEGsSapporo031DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo031DEGs,
                                                       universe = genelist_rnaseq,
                                                       OrgDb = organism, 
                                                       keyType = 'FLYBASE',
                                                       readable = T,
                                                       ont = "BP",
                                                       pvalueCutoff = 0.05, 
                                                       qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiNonDEGsSapporo031DEGs <- go_enrich_BP_SendaiNonDEGsSapporo031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiNonDEGsSapporo031DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/GO_BP_SendaiNonDEGsSapporo031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiNonDEGsSapporo031DEGs <- upsetplot(go_enrich_BP_SendaiNonDEGsSapporo031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/UpsetPlot_BP.pdf", g1_BP_SendaiNonDEGsSapporo031DEGs, h=6, w=10)

g2_BP_SendaiNonDEGsSapporo031DEGs <- barplot(go_enrich_BP_SendaiNonDEGsSapporo031DEGs, 
                                               drop = TRUE, 
                                               showCategory = 10, 
                                               title = "GO Biological Pathways",
                                               font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/BarPlot_BP.pdf", g2_BP_SendaiNonDEGsSapporo031DEGs, h=3, w=5)

g3_BP_SendaiNonDEGsSapporo031DEGs <- dotplot(go_enrich_BP_SendaiNonDEGsSapporo031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/DotPlot_BP.pdf", g3_BP_SendaiNonDEGsSapporo031DEGs, h=4, w=6)

g4_BP_SendaiNonDEGsSapporo031DEGs <- emapplot(go_enrich_BP_SendaiNonDEGsSapporo031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiNonDEGsSapporo031DEGs, h=6, w=8)

g5_BP_SendaiNonDEGsSapporo031DEGs <- goplot(go_enrich_BP_SendaiNonDEGsSapporo031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiNonDEGsSapporo031DEGs, h=6, w=8)

g6_BP_SendaiNonDEGsSapporo031DEGs <- cnetplot(go_enrich_BP_SendaiNonDEGsSapporo031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/CategoryNet_BP.pdf", g6_BP_SendaiNonDEGsSapporo031DEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiNonDEGsSapporo031DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo031DEGs,
                                                       universe = genelist_rnaseq,
                                                       OrgDb = organism, 
                                                       keyType = 'FLYBASE',
                                                       readable = T,
                                                       ont = "CC",
                                                       pvalueCutoff = 0.05, 
                                                       qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiNonDEGsSapporo031DEGs <- go_enrich_CC_SendaiNonDEGsSapporo031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiNonDEGsSapporo031DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/GO_CC_SendaiNonDEGsSapporo031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiNonDEGsSapporo031DEGs <- upsetplot(go_enrich_CC_SendaiNonDEGsSapporo031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/UpsetPlot_CC.pdf", g1_CC_SendaiNonDEGsSapporo031DEGs, h=6, w=10)

g2_CC_SendaiNonDEGsSapporo031DEGs <- barplot(go_enrich_CC_SendaiNonDEGsSapporo031DEGs, 
                                               drop = TRUE, 
                                               showCategory = 10, 
                                               title = "GO Cellular Components",
                                               font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/BarPlot_CC.pdf", g2_CC_SendaiNonDEGsSapporo031DEGs, h=3, w=5)

g3_CC_SendaiNonDEGsSapporo031DEGs <- dotplot(go_enrich_CC_SendaiNonDEGsSapporo031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/DotPlot_CC.pdf", g3_CC_SendaiNonDEGsSapporo031DEGs, h=4, w=5)

g4_CC_SendaiNonDEGsSapporo031DEGs <- emapplot(go_enrich_CC_SendaiNonDEGsSapporo031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiNonDEGsSapporo031DEGs, h=6, w=8)

g5_CC_SendaiNonDEGsSapporo031DEGs <- goplot(go_enrich_CC_SendaiNonDEGsSapporo031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiNonDEGsSapporo031DEGs, h=6, w=8)

g6_CC_SendaiNonDEGsSapporo031DEGs <- cnetplot(go_enrich_CC_SendaiNonDEGsSapporo031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/CategoryNet_CC.pdf", g6_CC_SendaiNonDEGsSapporo031DEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiNonDEGsSapporo031DEGs <- enrichGO(gene = genelist_SendaiNonDEGsSapporo031DEGs,
                                                       universe = genelist_rnaseq,
                                                       OrgDb = organism, 
                                                       keyType = 'FLYBASE',
                                                       readable = T,
                                                       ont = "MF",
                                                       pvalueCutoff = 0.05, 
                                                       qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiNonDEGsSapporo031DEGs <- go_enrich_MF_SendaiNonDEGsSapporo031DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiNonDEGsSapporo031DEGs, "../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/GO_MF_SendaiNonDEGsSapporo031DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiNonDEGsSapporo031DEGs <- upsetplot(go_enrich_MF_SendaiNonDEGsSapporo031DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/UpsetPlot_MF.pdf", g1_MF_SendaiNonDEGsSapporo031DEGs, h=6, w=10)

g2_MF_SendaiNonDEGsSapporo031DEGs <- barplot(go_enrich_MF_SendaiNonDEGsSapporo031DEGs, 
                                               drop = TRUE, 
                                               showCategory = 10, 
                                               title = "GO Molecular Functions",
                                               font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/BarPlot_MF.pdf", g2_MF_SendaiNonDEGsSapporo031DEGs, h=3, w=5)

g3_MF_SendaiNonDEGsSapporo031DEGs <- dotplot(go_enrich_MF_SendaiNonDEGsSapporo031DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/DotPlot_MF.pdf", g3_MF_SendaiNonDEGsSapporo031DEGs, h=4, w=6)

g4_MF_SendaiNonDEGsSapporo031DEGs <- emapplot(go_enrich_MF_SendaiNonDEGsSapporo031DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiNonDEGsSapporo031DEGs, h=6, w=8)

g5_MF_SendaiNonDEGsSapporo031DEGs <- goplot(go_enrich_MF_SendaiNonDEGsSapporo031DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiNonDEGsSapporo031DEGs, h=6, w=8)

g6_MF_SendaiNonDEGsSapporo031DEGs <- cnetplot(go_enrich_MF_SendaiNonDEGsSapporo031DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiNonDEGsSapporo031DEGs/CategoryNet_MF.pdf", g6_MF_SendaiNonDEGsSapporo031DEGs, h=12, w=10)



#### SendaiDEGsSapporo07NonDEGsEnrichedModule ####
genelist_SendaiDEGsSapporo07NonDEGsEnrichedModule = read.delim("../figures/tmp/TableS5_SendaiDEGsSapporo07NonDEGsEnrichedModule.tsv", header = T) 
genelist_SendaiDEGsSapporo07NonDEGsEnrichedModule <- as.vector(genelist_SendaiDEGsSapporo07NonDEGsEnrichedModule$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGsEnrichedModule,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "BP",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/GO_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- upsetplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/UpsetPlot_BP.pdf", g1_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=10)

g2_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- barplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Biological Pathways",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/BarPlot_BP.pdf", g2_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=3, w=5)

g3_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- dotplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/DotPlot_BP.pdf", g3_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=4, w=6)

g4_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- emapplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/EncrichmentMap_BP.pdf", g4_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=8)

g5_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- goplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/EncrichedGOGraph_BP.pdf", g5_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=8)

g6_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule <- cnetplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/CategoryNet_BP.pdf", g6_BP_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGsEnrichedModule,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "CC",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/GO_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- upsetplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/UpsetPlot_CC.pdf", g1_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=10)

g2_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- barplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Cellular Components",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/BarPlot_CC.pdf", g2_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=3, w=5)

g3_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- dotplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/DotPlot_CC.pdf", g3_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=4, w=5)

g4_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- emapplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/EncrichmentMap_CC.pdf", g4_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=8)

g5_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- goplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/EncrichedGOGraph_CC.pdf", g5_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=8)

g6_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule <- cnetplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/CategoryNet_CC.pdf", g6_CC_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGsEnrichedModule,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "MF",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/GO_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- upsetplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/UpsetPlot_MF.pdf", g1_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=10)

g2_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- barplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Molecular Functions",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/BarPlot_MF.pdf", g2_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=3, w=5)

g3_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- dotplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/DotPlot_MF.pdf", g3_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=4, w=6)

g4_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- emapplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/EncrichmentMap_MF.pdf", g4_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=8)

g5_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- goplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/EncrichedGOGraph_MF.pdf", g5_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=6, w=8)

g6_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule <- cnetplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGsEnrichedModule/CategoryNet_MF.pdf", g6_MF_SendaiDEGsSapporo07NonDEGsEnrichedModule, h=12, w=10)




#### SendaiDEGsSapporo07NonDEGs31DEGs ####
genelist_SendaiDEGsSapporo07NonDEGs31DEGs = read.delim("../figures/tmp/TableS5_SendaiDEGsSapporo07NonDEGs31DEGs.tsv", header = T) 
genelist_SendaiDEGsSapporo07NonDEGs31DEGs <- as.vector(genelist_SendaiDEGsSapporo07NonDEGs31DEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGs31DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "BP",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/GO_BP_SendaiDEGsSapporo07NonDEGs31DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- upsetplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/UpsetPlot_BP.pdf", g1_BP_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=10)

g2_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- barplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Biological Pathways",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/BarPlot_BP.pdf", g2_BP_SendaiDEGsSapporo07NonDEGs31DEGs, h=3, w=5)

g3_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- dotplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/DotPlot_BP.pdf", g3_BP_SendaiDEGsSapporo07NonDEGs31DEGs, h=4, w=6)

g4_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- emapplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=8)

g5_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- goplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=8)

g6_BP_SendaiDEGsSapporo07NonDEGs31DEGs <- cnetplot(go_enrich_BP_SendaiDEGsSapporo07NonDEGs31DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/CategoryNet_BP.pdf", g6_BP_SendaiDEGsSapporo07NonDEGs31DEGs, h=12, w=10)

##### CC #####
go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGs31DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "CC",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/GO_CC_SendaiDEGsSapporo07NonDEGs31DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- upsetplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/UpsetPlot_CC.pdf", g1_CC_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=10)

g2_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- barplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Cellular Components",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/BarPlot_CC.pdf", g2_CC_SendaiDEGsSapporo07NonDEGs31DEGs, h=3, w=5)

g3_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- dotplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/DotPlot_CC.pdf", g3_CC_SendaiDEGsSapporo07NonDEGs31DEGs, h=4, w=5)

g4_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- emapplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=8)

g5_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- goplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=8)

g6_CC_SendaiDEGsSapporo07NonDEGs31DEGs <- cnetplot(go_enrich_CC_SendaiDEGsSapporo07NonDEGs31DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/CategoryNet_CC.pdf", g6_CC_SendaiDEGsSapporo07NonDEGs31DEGs, h=12, w=10)

##### MF #####
go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- enrichGO(gene = genelist_SendaiDEGsSapporo07NonDEGs31DEGs,
                                                    universe = genelist_rnaseq,
                                                    OrgDb = organism, 
                                                    keyType = 'FLYBASE',
                                                    readable = T,
                                                    ont = "MF",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs, "../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/GO_MF_SendaiDEGsSapporo07NonDEGs31DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- upsetplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/UpsetPlot_MF.pdf", g1_MF_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=10)

g2_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- barplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs, 
                                            drop = TRUE, 
                                            showCategory = 10, 
                                            title = "GO Molecular Functions",
                                            font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/BarPlot_MF.pdf", g2_MF_SendaiDEGsSapporo07NonDEGs31DEGs, h=3, w=5)

g3_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- dotplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/DotPlot_MF.pdf", g3_MF_SendaiDEGsSapporo07NonDEGs31DEGs, h=4, w=6)

g4_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- emapplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=8)

g5_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- goplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiDEGsSapporo07NonDEGs31DEGs, h=6, w=8)

g6_MF_SendaiDEGsSapporo07NonDEGs31DEGs <- cnetplot(go_enrich_MF_SendaiDEGsSapporo07NonDEGs31DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/CategoryNet_MF.pdf", g6_MF_SendaiDEGsSapporo07NonDEGs31DEGs, h=12, w=10)



#### PopSpe07DEGs ####
genelist_PopSpe07DEGs = read.delim("../figures/tmp/TableS5_PopSpe07DEGs.tsv", header = T) 
genelist_PopSpe07DEGs <- as.vector(genelist_PopSpe07DEGs$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_PopSpe07DEGs <- enrichGO(gene = genelist_PopSpe07DEGs,
                                                          universe = genelist_rnaseq,
                                                          OrgDb = organism, 
                                                          keyType = 'FLYBASE',
                                                          readable = T,
                                                          ont = "BP",
                                                          pvalueCutoff = 0.05, 
                                                          qvalueCutoff = 0.10)

df_go_enrich_BP_PopSpe07DEGs <- go_enrich_BP_PopSpe07DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_PopSpe07DEGs, "../data/enrichment/clusterProfiler_PopSpe07DEGs/GO_BP_PopSpe07DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_PopSpe07DEGs <- upsetplot(go_enrich_BP_PopSpe07DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/UpsetPlot_BP.pdf", g1_BP_PopSpe07DEGs, h=6, w=10)

g2_BP_PopSpe07DEGs <- barplot(go_enrich_BP_PopSpe07DEGs, 
                                                  drop = TRUE, 
                                                  showCategory = 10, 
                                                  title = "GO Biological Pathways",
                                                  font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/BarPlot_BP.pdf", g2_BP_PopSpe07DEGs, h=3, w=5)

g3_BP_PopSpe07DEGs <- dotplot(go_enrich_BP_PopSpe07DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/DotPlot_BP.pdf", g3_BP_PopSpe07DEGs, h=4, w=6)

g4_BP_PopSpe07DEGs <- emapplot(go_enrich_BP_PopSpe07DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/EncrichmentMap_BP.pdf", g4_BP_PopSpe07DEGs, h=6, w=8)

g5_BP_PopSpe07DEGs <- goplot(go_enrich_BP_PopSpe07DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/EncrichedGOGraph_BP.pdf", g5_BP_PopSpe07DEGs, h=6, w=8)

g6_BP_PopSpe07DEGs <- cnetplot(go_enrich_BP_PopSpe07DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/CategoryNet_BP.pdf", g6_BP_PopSpe07DEGs, h=12, w=10)

##### CC #####
go_enrich_CC_PopSpe07DEGs <- enrichGO(gene = genelist_PopSpe07DEGs,
                                                          universe = genelist_rnaseq,
                                                          OrgDb = organism, 
                                                          keyType = 'FLYBASE',
                                                          readable = T,
                                                          ont = "CC",
                                                          pvalueCutoff = 0.05, 
                                                          qvalueCutoff = 0.10)

df_go_enrich_CC_PopSpe07DEGs <- go_enrich_CC_PopSpe07DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_PopSpe07DEGs, "../data/enrichment/clusterProfiler_PopSpe07DEGs/GO_CC_PopSpe07DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_PopSpe07DEGs <- upsetplot(go_enrich_CC_PopSpe07DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/UpsetPlot_CC.pdf", g1_CC_PopSpe07DEGs, h=6, w=10)

g2_CC_PopSpe07DEGs <- barplot(go_enrich_CC_PopSpe07DEGs, 
                                                  drop = TRUE, 
                                                  showCategory = 10, 
                                                  title = "GO Cellular Components",
                                                  font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/BarPlot_CC.pdf", g2_CC_PopSpe07DEGs, h=3, w=5)

g3_CC_PopSpe07DEGs <- dotplot(go_enrich_CC_PopSpe07DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/DotPlot_CC.pdf", g3_CC_PopSpe07DEGs, h=4, w=5)

g4_CC_PopSpe07DEGs <- emapplot(go_enrich_CC_PopSpe07DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/EncrichmentMap_CC.pdf", g4_CC_PopSpe07DEGs, h=6, w=8)

g5_CC_PopSpe07DEGs <- goplot(go_enrich_CC_PopSpe07DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/EncrichedGOGraph_CC.pdf", g5_CC_PopSpe07DEGs, h=6, w=8)

g6_CC_PopSpe07DEGs <- cnetplot(go_enrich_CC_PopSpe07DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/CategoryNet_CC.pdf", g6_CC_PopSpe07DEGs, h=12, w=10)

##### MF #####
go_enrich_MF_PopSpe07DEGs <- enrichGO(gene = genelist_PopSpe07DEGs,
                                                          universe = genelist_rnaseq,
                                                          OrgDb = organism, 
                                                          keyType = 'FLYBASE',
                                                          readable = T,
                                                          ont = "MF",
                                                          pvalueCutoff = 0.05, 
                                                          qvalueCutoff = 0.10)

df_go_enrich_MF_PopSpe07DEGs <- go_enrich_MF_PopSpe07DEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_PopSpe07DEGs, "../data/enrichment/clusterProfiler_PopSpe07DEGs/GO_MF_PopSpe07DEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_PopSpe07DEGs <- upsetplot(go_enrich_MF_PopSpe07DEGs)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/UpsetPlot_MF.pdf", g1_MF_PopSpe07DEGs, h=6, w=10)

g2_MF_PopSpe07DEGs <- barplot(go_enrich_MF_PopSpe07DEGs, 
                                                  drop = TRUE, 
                                                  showCategory = 10, 
                                                  title = "GO Molecular Functions",
                                                  font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/BarPlot_MF.pdf", g2_MF_PopSpe07DEGs, h=3, w=5)

g3_MF_PopSpe07DEGs <- dotplot(go_enrich_MF_PopSpe07DEGs, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/DotPlot_MF.pdf", g3_MF_PopSpe07DEGs, h=4, w=6)

g4_MF_PopSpe07DEGs <- emapplot(go_enrich_MF_PopSpe07DEGs %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/EncrichmentMap_MF.pdf", g4_MF_PopSpe07DEGs, h=6, w=8)

g5_MF_PopSpe07DEGs <- goplot(go_enrich_MF_PopSpe07DEGs, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/EncrichedGOGraph_MF.pdf", g5_MF_PopSpe07DEGs, h=6, w=8)

g6_MF_PopSpe07DEGs <- cnetplot(go_enrich_MF_PopSpe07DEGs, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_PopSpe07DEGs/CategoryNet_MF.pdf", g6_MF_PopSpe07DEGs, h=12, w=10)


