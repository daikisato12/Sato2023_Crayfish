#### load libraries ####
# rm(list = ls(all.names = TRUE))
# library(tidyverse)
# install.packages("BiocManager")
library(BiocManager)
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)

organism = "org.Dm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#### load genes ####
setwd("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/codes/")

# reading dataset
genelist_background = read.delim("../data/blastp/Allgenes_LG_crayfish_newID_flyID.txt", h=F)
genelist_background <- as.vector(genelist_background[,1]) %>% unique()

# reading dataset
genelist_rnaseq = read.delim("../data/RNA-seq/rawdata/RNAseq_analyzed_CrayfishID_flyID.txt", h=F)
genelist_rnaseq <- as.vector(genelist_rnaseq[,1]) %>% unique()

genelist_DEG2 = read.delim("../data/RNA-seq/DEGs/RNAseq_DEG2_CrayfishID_flyID.txt", header=F)
genelist_DEG2 <- as.vector(genelist_DEG2[,1]) %>% unique()

genelist_PBS = read.delim("../data/PBS/1kbp/PBS_mincov4_maxcov16_1kbp-top1_2combsharedloci_sharedgene_flyID.txt", header=F)
genelist_PBS <- as.vector(genelist_PBS[,1]) %>% unique()

genelist_PBS_DEG2_shared = read.delim("../data/PBS/1kbp/PBS_mincov4_maxcov16_1kbp-top1_2combsharedloci_DEGs_sharedgene_flyID.txt", header=F)
genelist_PBS_DEG2_shared <- as.vector(genelist_PBS_DEG2_shared[,1]) %>% unique()

# BiocManager::install("enrichplot", force = TRUE)
# install.packages("ggupset")
# install.packages("ggnewscale")
library(ggupset)
library(ggnewscale)
library(enrichplot)

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
genelist_Module = read.delim("../figures/tmp/WGCNA_Module.tsv", header = T) 
genelist_Module <- as.vector(genelist_Module$FlyID) %>% na.omit() %>% unique()

##### BP #####
go_enrich_BP_Module <- enrichGO(gene = genelist_Module,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "BP",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_BP_Module <- go_enrich_BP_Module@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_Module, "../data/enrichment/clusterProfiler_Module/GO_BP_Module.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_Module <- upsetplot(go_enrich_BP_Module)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/UpsetPlot_BP.pdf", g1_BP_Module, h=6, w=10)

g2_BP_Module <- barplot(go_enrich_BP_Module, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Biological Pathways",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/BarPlot_BP.pdf", g2_BP_Module, h=3, w=5)

g3_BP_Module <- dotplot(go_enrich_BP_Module, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/DotPlot_BP.pdf", g3_BP_Module, h=4, w=6)

g4_BP_Module <- emapplot(go_enrich_BP_Module %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/EncrichmentMap_BP.pdf", g4_BP_Module, h=6, w=8)

g5_BP_Module <- goplot(go_enrich_BP_Module, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/EncrichedGOGraph_BP.pdf", g5_BP_Module, h=6, w=8)

g6_BP_Module <- cnetplot(go_enrich_BP_Module, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/CategoryNet_BP.pdf", g6_BP_Module, h=12, w=10)

##### CC #####
go_enrich_CC_Module <- enrichGO(gene = genelist_Module,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "CC",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_CC_Module <- go_enrich_CC_Module@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_Module, "../data/enrichment/clusterProfiler_Module/GO_CC_Module.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_Module <- upsetplot(go_enrich_CC_Module)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/UpsetPlot_CC.pdf", g1_CC_Module, h=6, w=10)

g2_CC_Module <- barplot(go_enrich_CC_Module, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Cellular Components",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/BarPlot_CC.pdf", g2_CC_Module, h=3, w=5)

g3_CC_Module <- dotplot(go_enrich_CC_Module, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/DotPlot_CC.pdf", g3_CC_Module, h=4, w=5)

g4_CC_Module <- emapplot(go_enrich_CC_Module %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/EncrichmentMap_CC.pdf", g4_CC_Module, h=6, w=8)

g5_CC_Module <- goplot(go_enrich_CC_Module, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/EncrichedGOGraph_CC.pdf", g5_CC_Module, h=6, w=8)

g6_CC_Module <- cnetplot(go_enrich_CC_Module, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/CategoryNet_CC.pdf", g6_CC_Module, h=12, w=10)

##### MF #####
go_enrich_MF_Module <- enrichGO(gene = genelist_Module,
                              universe = genelist_rnaseq,
                              OrgDb = organism, 
                              keyType = 'FLYBASE',
                              readable = T,
                              ont = "MF",
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)

df_go_enrich_MF_Module <- go_enrich_MF_Module@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_Module, "../data/enrichment/clusterProfiler_Module/GO_MF_Module.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_Module <- upsetplot(go_enrich_MF_Module)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/UpsetPlot_MF.pdf", g1_MF_Module, h=6, w=10)

g2_MF_Module <- barplot(go_enrich_MF_Module, 
                      drop = TRUE, 
                      showCategory = 10, 
                      title = "GO Molecular Functions",
                      font.size = 8)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/BarPlot_MF.pdf", g2_MF_Module, h=3, w=5)

g3_MF_Module <- dotplot(go_enrich_MF_Module, orderBy = "x")
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/DotPlot_MF.pdf", g3_MF_Module, h=4, w=6)

g4_MF_Module <- emapplot(go_enrich_MF_Module %>% pairwise_termsim())
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/EncrichmentMap_MF.pdf", g4_MF_Module, h=6, w=8)

g5_MF_Module <- goplot(go_enrich_MF_Module, showCategory = 10)
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/EncrichedGOGraph_MF.pdf", g5_MF_Module, h=6, w=8)

g6_MF_Module <- cnetplot(go_enrich_MF_Module, categorySize="pvalue")
ggplot2::ggsave("../data/enrichment/clusterProfiler_Module/CategoryNet_MF.pdf", g6_MF_Module, h=12, w=10)

