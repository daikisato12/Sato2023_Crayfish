# rm(list = ls(all.names = TRUE))
# install.packages("tidyverse")
library(tidyverse)
# install.packages("BiocManager")
library(BiocManager)
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# install.packages("wordcloud")
library(clusterProfiler)
# library(wordcloud)

organism = "org.Dm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# BiocManager::install("enrichplot", force = TRUE)
# install.packages("ggupset")
# install.packages("ggnewscale")
library(ggupset)
# library(ggnewscale)
library(enrichplot)

# BiocManager::install("ReactomePA")
library(ReactomePA)
BiocManager::install("biomaRt")
library(biomaRt)

#### load genes ####
# loading genes analyzed in RNA-seq
genelist_rnaseq = read.delim("../../data/analyzed_data/RNAseq/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104_flyID.txt", h=F)
genelist_rnaseq <- as.vector(genelist_rnaseq[,1]) %>% unique()
genelist_rnaseq_entrez <- mapIds(org.Dm.eg.db, 
                                 keys=genelist_rnaseq, 
                                 column="ENTREZID", 
                                 keytype="ENSEMBL",
                                 multiVals="first")

mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))

#### analysis ####
##### Sendai-specific DEGs #####
dir.create("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs", recursive = TRUE)

genelist_SendaiSpecificDEGs = read.delim("../../data/analyzed_data/RNAseq/genes/DEG_SendaiSpecific.tsv", header = T) 
genelist_SendaiSpecificDEGs <- as.vector(genelist_SendaiSpecificDEGs$FlyID) %>% na.omit() %>% unique()
genelist_SendaiSpecificDEGs_entrez <- mapIds(org.Dm.eg.db, 
                                         keys=genelist_SendaiSpecificDEGs, 
                                         column="ENTREZID", 
                                         keytype="ENSEMBL",
                                         multiVals="first")

###### KEGG ######
KEGG_enrich_SendaiSpecificDEGs <- enrichKEGG(gene = genelist_SendaiSpecificDEGs_entrez, 
                                         universe = genelist_rnaseq_entrez,
                                         organism = "dme", 
                                         pvalueCutoff = 0.05,
                                         keyType = "ncbi-geneid",
                                         # minGSSize = 10, 
                                         # maxGSSize = 500, 
                                         # pAdjustMethod = "BH", 
                                         qvalueCutoff = 0.10)

genelist_SendaiSpecificDEGs_entrez_names <- genelist_SendaiSpecificDEGs_entrez
genelist_SendaiSpecificDEGs_entrez_flyid <- names(genelist_SendaiSpecificDEGs_entrez)
names(genelist_SendaiSpecificDEGs_entrez_flyid) <- genelist_SendaiSpecificDEGs_entrez_names
df_KEGG_enrich_SendaiSpecificDEGs <- KEGG_enrich_SendaiSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", ")) %>%
  separate_rows(geneID, sep = ", ") %>%
  mutate(geneID = genelist_SendaiSpecificDEGs_entrez_flyid[geneID])

df_KEGG_enrich_SendaiSpecificDEGs_genename <- getBM(filters = "ensembl_gene_id", 
                                                     attributes = c("ensembl_gene_id", "external_gene_name"),
                                                     values = df_KEGG_enrich_SendaiSpecificDEGs$geneID, mart = mart) %>%
  rename(geneID = ensembl_gene_id,
         genename = external_gene_name)

df_KEGG_enrich_SendaiSpecificDEGs2 <- left_join(df_KEGG_enrich_SendaiSpecificDEGs,
                                                 df_KEGG_enrich_SendaiSpecificDEGs_genename) %>%
group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
  summarize(geneID = str_c(geneID, collapse = ", "),
            genename = str_c(genename, collapse = ", "),)
write.table(df_KEGG_enrich_SendaiSpecificDEGs2, "../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/KEGG_enrichment_SendaiSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


KEGG_enrich_SendaiSpecificDEGs@result
barplot(KEGG_enrich_SendaiSpecificDEGs, showCategory=10)
browseKEGG(KEGG_enrich_SendaiSpecificDEGs, "dme00514")


###### Reactome ######
Reactome_enrich_SendaiSpecificDEGs <- enrichPathway(gene = genelist_SendaiSpecificDEGs_entrez, 
                                                universe = genelist_rnaseq_entrez,
                                                organism = "fly",
                                                pvalueCutoff = 0.05, 
                                                qvalueCutoff = 0.1,
                                                readable = T)

df_Reactome_enrich_SendaiSpecificDEGs <- Reactome_enrich_SendaiSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_Reactome_enrich_SendaiSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/Reactome_enrichment_SendaiSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


###### GO BP ######
go_enrich_BP_SendaiSpecificDEGs <- enrichGO(gene = genelist_SendaiSpecificDEGs,
                                        universe = genelist_rnaseq,
                                        OrgDb = organism, 
                                        keyType = 'FLYBASE',
                                        readable = T,
                                        ont = "BP",
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 0.10)

df_go_enrich_BP_SendaiSpecificDEGs <- go_enrich_BP_SendaiSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SendaiSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/GO_BP_SendaiSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SendaiSpecificDEGs <- upsetplot(go_enrich_BP_SendaiSpecificDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/UpsetPlot_BP.pdf", g1_BP_SendaiSpecificDEGs, h=6, w=10)

g2_BP_SendaiSpecificDEGs <- barplot(go_enrich_BP_SendaiSpecificDEGs, 
                                drop = TRUE, 
                                showCategory = 10, 
                                title = "GO Biological Pathways",
                                font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/BarPlot_BP.pdf", g2_BP_SendaiSpecificDEGs, h=3, w=5)

g3_BP_SendaiSpecificDEGs <- dotplot(go_enrich_BP_SendaiSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/DotPlot_BP.pdf", g3_BP_SendaiSpecificDEGs, h=4, w=6)

g4_BP_SendaiSpecificDEGs <- emapplot(go_enrich_BP_SendaiSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/EncrichmentMap_BP.pdf", g4_BP_SendaiSpecificDEGs, h=6, w=8)

g5_BP_SendaiSpecificDEGs <- goplot(go_enrich_BP_SendaiSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SendaiSpecificDEGs, h=6, w=8)

g6_BP_SendaiSpecificDEGs <- cnetplot(go_enrich_BP_SendaiSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/CategoryNet_BP.pdf", g6_BP_SendaiSpecificDEGs, h=12, w=10)


###### GO CC ######
go_enrich_CC_SendaiSpecificDEGs <- enrichGO(gene = genelist_SendaiSpecificDEGs,
                                        universe = genelist_rnaseq,
                                        OrgDb = organism, 
                                        keyType = 'FLYBASE',
                                        readable = T,
                                        ont = "CC",
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 0.10)

df_go_enrich_CC_SendaiSpecificDEGs <- go_enrich_CC_SendaiSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SendaiSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/GO_CC_SendaiSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SendaiSpecificDEGs <- upsetplot(go_enrich_CC_SendaiSpecificDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/UpsetPlot_CC.pdf", g1_CC_SendaiSpecificDEGs, h=6, w=10)

g2_CC_SendaiSpecificDEGs <- barplot(go_enrich_CC_SendaiSpecificDEGs, 
                                drop = TRUE, 
                                showCategory = 10, 
                                title = "GO Cellular Components",
                                font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/BarPlot_CC.pdf", g2_CC_SendaiSpecificDEGs, h=3, w=5)

g3_CC_SendaiSpecificDEGs <- dotplot(go_enrich_CC_SendaiSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/DotPlot_CC.pdf", g3_CC_SendaiSpecificDEGs, h=4, w=5)

g4_CC_SendaiSpecificDEGs <- emapplot(go_enrich_CC_SendaiSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/EncrichmentMap_CC.pdf", g4_CC_SendaiSpecificDEGs, h=6, w=8)

g5_CC_SendaiSpecificDEGs <- goplot(go_enrich_CC_SendaiSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SendaiSpecificDEGs, h=6, w=8)

g6_CC_SendaiSpecificDEGs <- cnetplot(go_enrich_CC_SendaiSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/CategoryNet_CC.pdf", g6_CC_SendaiSpecificDEGs, h=12, w=10)

###### GO MF ######
go_enrich_MF_SendaiSpecificDEGs <- enrichGO(gene = genelist_SendaiSpecificDEGs,
                                        universe = genelist_rnaseq,
                                        OrgDb = organism, 
                                        keyType = 'FLYBASE',
                                        readable = T,
                                        ont = "MF",
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 0.10)

df_go_enrich_MF_SendaiSpecificDEGs <- go_enrich_MF_SendaiSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SendaiSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/GO_MF_SendaiSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SendaiSpecificDEGs <- upsetplot(go_enrich_MF_SendaiSpecificDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/UpsetPlot_MF.pdf", g1_MF_SendaiSpecificDEGs, h=6, w=10)

g2_MF_SendaiSpecificDEGs <- barplot(go_enrich_MF_SendaiSpecificDEGs, 
                                drop = TRUE, 
                                showCategory = 10, 
                                title = "GO Molecular Functions",
                                font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/BarPlot_MF.pdf", g2_MF_SendaiSpecificDEGs, h=3, w=5)

g3_MF_SendaiSpecificDEGs <- dotplot(go_enrich_MF_SendaiSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/DotPlot_MF.pdf", g3_MF_SendaiSpecificDEGs, h=4, w=6)

g4_MF_SendaiSpecificDEGs <- emapplot(go_enrich_MF_SendaiSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/EncrichmentMap_MF.pdf", g4_MF_SendaiSpecificDEGs, h=6, w=8)

g5_MF_SendaiSpecificDEGs <- goplot(go_enrich_MF_SendaiSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SendaiSpecificDEGs, h=6, w=8)

g6_MF_SendaiSpecificDEGs <- cnetplot(go_enrich_MF_SendaiSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/CategoryNet_MF.pdf", g6_MF_SendaiSpecificDEGs, h=12, w=10)



##### Sapporo-specific DEGs #####
dir.create("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs", recursive = TRUE)

genelist_SapporoSpecificDEGs = read.delim("../../data/analyzed_data/RNAseq/genes/DEG_SapporoSpecific.tsv", header = T) 
genelist_SapporoSpecificDEGs <- as.vector(genelist_SapporoSpecificDEGs$FlyID) %>% na.omit() %>% unique()
genelist_SapporoSpecificDEGs_entrez <- mapIds(org.Dm.eg.db, 
                                             keys=genelist_SapporoSpecificDEGs, 
                                             column="ENTREZID", 
                                             keytype="ENSEMBL",
                                             multiVals="first")

###### KEGG ######
KEGG_enrich_SapporoSpecificDEGs <- enrichKEGG(gene = genelist_SapporoSpecificDEGs_entrez, 
                                             universe = genelist_rnaseq_entrez,
                                             organism = "dme", 
                                             pvalueCutoff = 0.05,
                                             keyType = "ncbi-geneid",
                                             # minGSSize = 10, 
                                             # maxGSSize = 500, 
                                             # pAdjustMethod = "BH", 
                                             qvalueCutoff = 0.10)

genelist_SapporoSpecificDEGs_entrez_names <- genelist_SapporoSpecificDEGs_entrez
genelist_SapporoSpecificDEGs_entrez_flyid <- names(genelist_SapporoSpecificDEGs_entrez)
names(genelist_SapporoSpecificDEGs_entrez_flyid) <- genelist_SapporoSpecificDEGs_entrez_names

df_KEGG_enrich_SapporoSpecificDEGs <- KEGG_enrich_SapporoSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", ")) %>%
  separate_rows(geneID, sep = ", ") %>%
  mutate(geneID = genelist_SapporoSpecificDEGs_entrez_flyid[geneID])

df_KEGG_enrich_SapporoSpecificDEGs_genename <- getBM(filters = "ensembl_gene_id", 
      attributes = c("ensembl_gene_id", "external_gene_name"),
      values = df_KEGG_enrich_SapporoSpecificDEGs$geneID, mart = mart) %>%
  rename(geneID = ensembl_gene_id,
         genename = external_gene_name)

df_KEGG_enrich_SapporoSpecificDEGs2 <- left_join(df_KEGG_enrich_SapporoSpecificDEGs,
                                                 df_KEGG_enrich_SapporoSpecificDEGs_genename) %>%
  group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
  summarize(geneID = str_c(geneID, collapse = ", "),
            genename = str_c(genename, collapse = ", "),)
write.table(df_KEGG_enrich_SapporoSpecificDEGs2, "../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/KEGG_enrichment_SapporoSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


###### Reactome ######
Reactome_enrich_SapporoSpecificDEGs <- enrichPathway(gene = genelist_SapporoSpecificDEGs_entrez, 
                                                    universe = genelist_rnaseq_entrez,
                                                    organism = "fly",
                                                    pvalueCutoff = 0.05, 
                                                    qvalueCutoff = 0.1,
                                                    readable = T)

df_Reactome_enrich_SapporoSpecificDEGs <- Reactome_enrich_SapporoSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_Reactome_enrich_SapporoSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/Reactome_enrichment_SapporoSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


###### GO BP ######
go_enrich_BP_SapporoSpecificDEGs <- enrichGO(gene = genelist_SapporoSpecificDEGs,
                                            universe = genelist_rnaseq,
                                            OrgDb = organism, 
                                            keyType = 'FLYBASE',
                                            readable = T,
                                            ont = "BP",
                                            pvalueCutoff = 0.05, 
                                            qvalueCutoff = 0.10)

df_go_enrich_BP_SapporoSpecificDEGs <- go_enrich_BP_SapporoSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_go_enrich_BP_SapporoSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/GO_BP_SapporoSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SapporoSpecificDEGs <- upsetplot(go_enrich_BP_SapporoSpecificDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/UpsetPlot_BP.pdf", g1_BP_SapporoSpecificDEGs, h=6, w=10)

g2_BP_SapporoSpecificDEGs <- barplot(go_enrich_BP_SapporoSpecificDEGs, 
                                    drop = TRUE, 
                                    showCategory = 10, 
                                    title = "GO Biological Pathways",
                                    font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/BarPlot_BP.pdf", g2_BP_SapporoSpecificDEGs, h=3, w=5)

g3_BP_SapporoSpecificDEGs <- dotplot(go_enrich_BP_SapporoSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/DotPlot_BP.pdf", g3_BP_SapporoSpecificDEGs, h=4, w=6)

g4_BP_SapporoSpecificDEGs <- emapplot(go_enrich_BP_SapporoSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/EncrichmentMap_BP.pdf", g4_BP_SapporoSpecificDEGs, h=6, w=8)

g5_BP_SapporoSpecificDEGs <- goplot(go_enrich_BP_SapporoSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SapporoSpecificDEGs, h=6, w=8)

g6_BP_SapporoSpecificDEGs <- cnetplot(go_enrich_BP_SapporoSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/CategoryNet_BP.pdf", g6_BP_SapporoSpecificDEGs, h=12, w=10)


###### GO CC ######
go_enrich_CC_SapporoSpecificDEGs <- enrichGO(gene = genelist_SapporoSpecificDEGs,
                                            universe = genelist_rnaseq,
                                            OrgDb = organism, 
                                            keyType = 'FLYBASE',
                                            readable = T,
                                            ont = "CC",
                                            pvalueCutoff = 0.05, 
                                            qvalueCutoff = 0.10)

df_go_enrich_CC_SapporoSpecificDEGs <- go_enrich_CC_SapporoSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_go_enrich_CC_SapporoSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/GO_CC_SapporoSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SapporoSpecificDEGs <- upsetplot(go_enrich_CC_SapporoSpecificDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/UpsetPlot_CC.pdf", g1_CC_SapporoSpecificDEGs, h=6, w=10)

g2_CC_SapporoSpecificDEGs <- barplot(go_enrich_CC_SapporoSpecificDEGs, 
                                    drop = TRUE, 
                                    showCategory = 10, 
                                    title = "GO Cellular Components",
                                    font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/BarPlot_CC.pdf", g2_CC_SapporoSpecificDEGs, h=3, w=5)

g3_CC_SapporoSpecificDEGs <- dotplot(go_enrich_CC_SapporoSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/DotPlot_CC.pdf", g3_CC_SapporoSpecificDEGs, h=4, w=5)

g4_CC_SapporoSpecificDEGs <- emapplot(go_enrich_CC_SapporoSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/EncrichmentMap_CC.pdf", g4_CC_SapporoSpecificDEGs, h=6, w=8)

g5_CC_SapporoSpecificDEGs <- goplot(go_enrich_CC_SapporoSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SapporoSpecificDEGs, h=6, w=8)

g6_CC_SapporoSpecificDEGs <- cnetplot(go_enrich_CC_SapporoSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/CategoryNet_CC.pdf", g6_CC_SapporoSpecificDEGs, h=12, w=10)

###### GO MF ######
go_enrich_MF_SapporoSpecificDEGs <- enrichGO(gene = genelist_SapporoSpecificDEGs,
                                            universe = genelist_rnaseq,
                                            OrgDb = organism, 
                                            keyType = 'FLYBASE',
                                            readable = T,
                                            ont = "MF",
                                            pvalueCutoff = 0.05, 
                                            qvalueCutoff = 0.10)

df_go_enrich_MF_SapporoSpecificDEGs <- go_enrich_MF_SapporoSpecificDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_go_enrich_MF_SapporoSpecificDEGs, "../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/GO_MF_SapporoSpecificDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SapporoSpecificDEGs <- upsetplot(go_enrich_MF_SapporoSpecificDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/UpsetPlot_MF.pdf", g1_MF_SapporoSpecificDEGs, h=6, w=10)

g2_MF_SapporoSpecificDEGs <- barplot(go_enrich_MF_SapporoSpecificDEGs, 
                                    drop = TRUE, 
                                    showCategory = 10, 
                                    title = "GO Molecular Functions",
                                    font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/BarPlot_MF.pdf", g2_MF_SapporoSpecificDEGs, h=3, w=5)

g3_MF_SapporoSpecificDEGs <- dotplot(go_enrich_MF_SapporoSpecificDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/DotPlot_MF.pdf", g3_MF_SapporoSpecificDEGs, h=4, w=6)

g4_MF_SapporoSpecificDEGs <- emapplot(go_enrich_MF_SapporoSpecificDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/EncrichmentMap_MF.pdf", g4_MF_SapporoSpecificDEGs, h=6, w=8)

g5_MF_SapporoSpecificDEGs <- goplot(go_enrich_MF_SapporoSpecificDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SapporoSpecificDEGs, h=6, w=8)

g6_MF_SapporoSpecificDEGs <- cnetplot(go_enrich_MF_SapporoSpecificDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/CategoryNet_MF.pdf", g6_MF_SapporoSpecificDEGs, h=12, w=10)



##### Shared DEGs #####
dir.create("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs", recursive = TRUE)

genelist_SharedDEGs = read.delim("../../data/analyzed_data/RNAseq/genes/DEG_Shared.tsv", header = T) 
genelist_SharedDEGs <- as.vector(genelist_SharedDEGs$FlyID) %>% na.omit() %>% unique()
genelist_SharedDEGs_entrez <- mapIds(org.Dm.eg.db, 
                                              keys=genelist_SharedDEGs, 
                                              column="ENTREZID", 
                                              keytype="ENSEMBL",
                                              multiVals="first")

###### KEGG ######
KEGG_enrich_SharedDEGs <- enrichKEGG(gene = genelist_SharedDEGs_entrez, 
                                              universe = genelist_rnaseq_entrez,
                                              organism = "dme", 
                                              pvalueCutoff = 0.05,
                                              keyType = "ncbi-geneid",
                                              # minGSSize = 10, 
                                              # maxGSSize = 500, 
                                              # pAdjustMethod = "BH", 
                                              qvalueCutoff = 0.10)

genelist_SharedDEGs_entrez_names <- genelist_SharedDEGs_entrez
genelist_SharedDEGs_entrez_flyid <- names(genelist_SharedDEGs_entrez)
names(genelist_SharedDEGs_entrez_flyid) <- genelist_SharedDEGs_entrez_names
df_KEGG_enrich_SharedDEGs <- KEGG_enrich_SharedDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", ")) %>%
  separate_rows(geneID, sep = ", ") %>%
  mutate(geneID = genelist_SharedDEGs_entrez_flyid[geneID])

df_KEGG_enrich_SharedDEGs_genename <- getBM(filters = "ensembl_gene_id", 
                                                     attributes = c("ensembl_gene_id", "external_gene_name"),
                                                     values = df_KEGG_enrich_SharedDEGs$geneID, mart = mart) %>%
  rename(geneID = ensembl_gene_id,
         genename = external_gene_name)

df_KEGG_enrich_SharedDEGs2 <- left_join(df_KEGG_enrich_SharedDEGs,
                                                 df_KEGG_enrich_SharedDEGs_genename) %>%
group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
  summarize(geneID = str_c(geneID, collapse = ", "),
            genename = str_c(genename, collapse = ", "),)
write.table(df_KEGG_enrich_SharedDEGs2, "../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/KEGG_enrichment_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


###### Reactome ######
Reactome_enrich_SharedDEGs <- enrichPathway(gene = genelist_SharedDEGs_entrez, 
                                                     universe = genelist_rnaseq_entrez,
                                                     organism = "fly",
                                                     pvalueCutoff = 0.05, 
                                                     qvalueCutoff = 0.1,
                                                     readable = T)

df_Reactome_enrich_SharedDEGs <- Reactome_enrich_SharedDEGs@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
write.table(df_Reactome_enrich_SharedDEGs, "../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/Reactome_enrichment_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


###### GO BP ######
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
write.table(df_go_enrich_BP_SharedDEGs, "../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/GO_BP_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_SharedDEGs <- upsetplot(go_enrich_BP_SharedDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/UpsetPlot_BP.pdf", g1_BP_SharedDEGs, h=6, w=10)

g2_BP_SharedDEGs <- barplot(go_enrich_BP_SharedDEGs, 
                                     drop = TRUE, 
                                     showCategory = 10, 
                                     title = "GO Biological Pathways",
                                     font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/BarPlot_BP.pdf", g2_BP_SharedDEGs, h=3, w=5)

g3_BP_SharedDEGs <- dotplot(go_enrich_BP_SharedDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/DotPlot_BP.pdf", g3_BP_SharedDEGs, h=4, w=6)

g4_BP_SharedDEGs <- emapplot(go_enrich_BP_SharedDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/EncrichmentMap_BP.pdf", g4_BP_SharedDEGs, h=6, w=8)

g5_BP_SharedDEGs <- goplot(go_enrich_BP_SharedDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/EncrichedGOGraph_BP.pdf", g5_BP_SharedDEGs, h=6, w=8)

g6_BP_SharedDEGs <- cnetplot(go_enrich_BP_SharedDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/CategoryNet_BP.pdf", g6_BP_SharedDEGs, h=12, w=10)


###### GO CC ######
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
write.table(df_go_enrich_CC_SharedDEGs, "../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/GO_CC_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_SharedDEGs <- upsetplot(go_enrich_CC_SharedDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/UpsetPlot_CC.pdf", g1_CC_SharedDEGs, h=6, w=10)

g2_CC_SharedDEGs <- barplot(go_enrich_CC_SharedDEGs, 
                                     drop = TRUE, 
                                     showCategory = 10, 
                                     title = "GO Cellular Components",
                                     font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/BarPlot_CC.pdf", g2_CC_SharedDEGs, h=3, w=5)

g3_CC_SharedDEGs <- dotplot(go_enrich_CC_SharedDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/DotPlot_CC.pdf", g3_CC_SharedDEGs, h=4, w=5)

g4_CC_SharedDEGs <- emapplot(go_enrich_CC_SharedDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/EncrichmentMap_CC.pdf", g4_CC_SharedDEGs, h=6, w=8)

g5_CC_SharedDEGs <- goplot(go_enrich_CC_SharedDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/EncrichedGOGraph_CC.pdf", g5_CC_SharedDEGs, h=6, w=8)

g6_CC_SharedDEGs <- cnetplot(go_enrich_CC_SharedDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/CategoryNet_CC.pdf", g6_CC_SharedDEGs, h=12, w=10)

###### GO MF ######
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
write.table(df_go_enrich_MF_SharedDEGs, "../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/GO_MF_SharedDEGs.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_SharedDEGs <- upsetplot(go_enrich_MF_SharedDEGs)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/UpsetPlot_MF.pdf", g1_MF_SharedDEGs, h=6, w=10)

g2_MF_SharedDEGs <- barplot(go_enrich_MF_SharedDEGs, 
                                     drop = TRUE, 
                                     showCategory = 10, 
                                     title = "GO Molecular Functions",
                                     font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/BarPlot_MF.pdf", g2_MF_SharedDEGs, h=3, w=5)

g3_MF_SharedDEGs <- dotplot(go_enrich_MF_SharedDEGs, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/DotPlot_MF.pdf", g3_MF_SharedDEGs, h=4, w=6)

g4_MF_SharedDEGs <- emapplot(go_enrich_MF_SharedDEGs %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/EncrichmentMap_MF.pdf", g4_MF_SharedDEGs, h=6, w=8)

g5_MF_SharedDEGs <- goplot(go_enrich_MF_SharedDEGs, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/EncrichedGOGraph_MF.pdf", g5_MF_SharedDEGs, h=6, w=8)

g6_MF_SharedDEGs <- cnetplot(go_enrich_MF_SharedDEGs, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/CategoryNet_MF.pdf", g6_MF_SharedDEGs, h=12, w=10)

