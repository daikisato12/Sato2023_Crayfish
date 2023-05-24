#### load libraries ####
library(tidyverse)
# install.packages("ggrepel")
library(ggrepel)
library(ape)
# install.packages("phytools")
library(phytools)
# install.packages("geiger")
library(geiger)
library(nlme)

# install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# install.packages("wordcloud")
library(clusterProfiler)
library(wordcloud)
# install.packages("ggupset")
library(ggupset)
library(ggnewscale)
library(enrichplot)

organism = "org.Dm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

library(tidytext)


#### load dataset ####
df_prop <- read.delim("../../data/analyzed_data/PD/results_20230519.txt") %>%
  mutate(Group = if_else(Group == 1, "Non-invasive", "Invasive") %>%
           factor(., levels = c("Invasive", "Non-invasive")),
         Decapods = if_else(Decapods == 1, "Decapods", "Others") %>%
           factor(., levels = c("Decapods", "Others")),
         spe_abb = paste0(substr(Species, 1, 1), ". ", str_split(Species, " ") %>% map_chr(2))) %>%
  filter(!Species %in% c("Chinoecetes opilio", "Procambarus virginalis"))
colnames(df_prop) <- c("Invasive", "Decapods", "species", "phyname", "PD", "propsize", "spe_abb")

#### Figure 4a PD against propagule size ####
g4_1 <- ggplot(df_prop, aes(x = log10(propsize), y = PD)) +
  geom_smooth(method = 'lm', formula= y ~ x, color = "grey") +
  geom_point(aes(shape = Decapods, color = Invasive), size = 5) +
  ggrepel::geom_text_repel(aes(label = spe_abb),  fontface = "italic") +
  scale_color_manual(values = c("darkred", "grey")) +
  scale_shape_manual(values = c(17, 16)) +
  xlab(expression(Log[10]("Propagule size"))) +
  ylab(expression(italic(P)[D])) +
  theme_bw() +
  theme(legend.title = element_blank())

g4_1

#### Figure 4b residuals ####
##### all invertebrates #####
# fit model
model <- lm(PD ~ log10(propsize), data=df_prop)
# view model summary
# summary(model)

df_prop2 <- df_prop %>%
  bind_cols(data.frame(res_std = rstandard(model)))

g4_2_all <- ggplot(df_prop2, aes(x = res_std, y = reorder(species, res_std))) +
  geom_bar(stat = "identity", aes(fill = Invasive)) +
  # geom_text_repel(aes(label = species), fontface = "italic") + #nudge_x = 0.1, 
  xlab(expression(paste("Residuals of ", {italic(P)[D]}, " against propagule size", sep = " "))) +
  scale_fill_manual(values = c("darkred", "grey")) +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(face = "italic"),#element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())


###### PGLS analysis ######
# Phylogenetic tree
allTree <- ape::read.tree("../../data/analyzed_data/PD/phylogenetic_tree/withoutPvirCopi/iqtree_all.nwk")
phy <- multi2di(allTree) #make rooted tree
plot(phy)
# # Calculate PICs
df_all_pgls <- df_prop %>%
  tibble::column_to_rownames(var = "phyname")

pglsModel <- gls(PD ~ propsize, correlation = corBrownian(phy = allTree),
                 data = df_all_pgls, method = "ML")
summary(pglsModel)
plot(df_all_pgls$PD ~ df_all_pgls$propsize)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])


g4 <- ggpubr::ggarrange(g4_1, g4_2_all, ncol = 2, align = "h")
ggplot2::ggsave("../../figures/tmp/Figure4ab.pdf", g4, w = 14, h = 6)


#### Figure 4c number of rapidly evolving genes ####
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

#### Figure 4d enrichment analysis ####
list_cafe_all_genes <- read.delim("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/Gamma_change.tab", header = TRUE) %>%
  dplyr::rename(GeneID = FamilyID) %>%
  dplyr::select(!c(X.14.)) %>%
  pull(GeneID)
list_cafe_Pcla_expanded_genes <- df_cafe_Pcla_exp %>%
  pull(GeneID)

dir.create("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp", recursive = T)

##### BP #####
go_enrich_BP_Pcla_exp <- enrichGO(gene = list_cafe_Pcla_expanded_genes,
                                  universe = list_cafe_all_genes,
                                  OrgDb = organism, 
                                  keyType = 'FLYBASE',
                                  readable = T,
                                  ont = "BP",
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.10)

df_go_enrich_BP_Pcla_exp <- go_enrich_BP_Pcla_exp@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
if(nrow(df_go_enrich_BP_Pcla_exp) > 0){
  write.table(df_go_enrich_BP_Pcla_exp, "../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/GO_BP_Pcla_exp.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
}

g1_BP_Pcla_exp <- upsetplot(go_enrich_BP_Pcla_exp)
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/UpsetPlot_BP.pdf", g1_BP_Pcla_exp, h=6, w=10)

g2_BP_Pcla_exp <- barplot(go_enrich_BP_Pcla_exp, 
                          drop = TRUE, 
                          showCategory = 10, 
                          title = "GO Biological Pathways",
                          font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/BarPlot_BP.pdf", g2_BP_Pcla_exp, h=3, w=5)

g3_BP_Pcla_exp <- dotplot(go_enrich_BP_Pcla_exp, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/DotPlot_BP.pdf", g3_BP_Pcla_exp, h=4, w=6)

g4_BP_Pcla_exp <- emapplot(go_enrich_BP_Pcla_exp %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/EncrichmentMap_BP.pdf", g4_BP_Pcla_exp, h=6, w=8)

g4_BP_Pcla_exp <- goplot(go_enrich_BP_Pcla_exp, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/EncrichedGOGraph_BP.pdf", g4_BP_Pcla_exp, h=6, w=8)

g6_BP_Pcla_exp <- cnetplot(go_enrich_BP_Pcla_exp, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/CategoryNet_BP.pdf", g6_BP_Pcla_exp, h=12, w=10)

##### CC #####
# go_enrich_CC_Pcla_exp <- enrichGO(gene = list_cafe_Pcla_expanded_genes,
#                                   universe = list_cafe_all_genes,
#                                   OrgDb = organism, 
#                                   keyType = 'FLYBASE',
#                                   readable = T,
#                                   ont = "CC",
#                                   pvalueCutoff = 0.05, 
#                                   qvalueCutoff = 0.10)
# 
# df_go_enrich_CC_Pcla_exp <- go_enrich_CC_Pcla_exp@result %>%
#   filter(p.adjust < 0.05) %>%
#   mutate(geneID = str_replace_all(geneID, "/", ", "))
# if(nrow(df_go_enrich_CC_Pcla_exp) > 0){
#   write.table(df_go_enrich_CC_Pcla_exp, "../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/GO_CC_Pcla_exp.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
# }
# g1_CC_Pcla_exp <- upsetplot(go_enrich_CC_Pcla_exp)
# ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/UpsetPlot_CC.pdf", g1_CC_Pcla_exp, h=6, w=10)
# 
# g2_CC_Pcla_exp <- barplot(go_enrich_CC_Pcla_exp, 
#                           drop = TRUE, 
#                           showCategory = 10, 
#                           title = "GO Cellular Components",
#                           font.size = 8)
# ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/BarPlot_CC.pdf", g2_CC_Pcla_exp, h=3, w=5)
# 
# g3_CC_Pcla_exp <- dotplot(go_enrich_CC_Pcla_exp, orderBy = "x")
# ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/DotPlot_CC.pdf", g3_CC_Pcla_exp, h=4, w=5)
# 
# g4_CC_Pcla_exp <- emapplot(go_enrich_CC_Pcla_exp %>% pairwise_termsim())
# ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/EncrichmentMap_CC.pdf", g4_CC_Pcla_exp, h=6, w=8)
# 
# g4_CC_Pcla_exp <- goplot(go_enrich_CC_Pcla_exp, showCategory = 10)
# ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/EncrichedGOGraph_CC.pdf", g4_CC_Pcla_exp, h=6, w=8)
# 
# g6_CC_Pcla_exp <- cnetplot(go_enrich_CC_Pcla_exp, categorySize="pvalue")
# ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/CategoryNet_CC.pdf", g6_CC_Pcla_exp, h=12, w=10)

##### MF #####
go_enrich_MF_Pcla_exp <- enrichGO(gene = list_cafe_Pcla_expanded_genes,
                                  universe = list_cafe_all_genes,
                                  OrgDb = organism, 
                                  keyType = 'FLYBASE',
                                  readable = T,
                                  ont = "MF",
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.10)

df_go_enrich_MF_Pcla_exp <- go_enrich_MF_Pcla_exp@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "))
if(nrow(df_go_enrich_MF_Pcla_exp) > 0){
  write.table(df_go_enrich_MF_Pcla_exp, "../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/GO_MF_Pcla_exp.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
}

g1_MF_Pcla_exp <- upsetplot(go_enrich_MF_Pcla_exp)
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/UpsetPlot_MF.pdf", g1_MF_Pcla_exp, h=6, w=10)

g2_MF_Pcla_exp <- barplot(go_enrich_MF_Pcla_exp, 
                          drop = TRUE, 
                          showCategory = 10, 
                          title = "GO Molecular Functions",
                          font.size = 8)
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/BarPlot_MF.pdf", g2_MF_Pcla_exp, h=3, w=5)

g3_MF_Pcla_exp <- dotplot(go_enrich_MF_Pcla_exp, orderBy = "x")
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/DotPlot_MF.pdf", g3_MF_Pcla_exp, h=4, w=6)

g4_MF_Pcla_exp <- emapplot(go_enrich_MF_Pcla_exp %>% pairwise_termsim())
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/EncrichmentMap_MF.pdf", g4_MF_Pcla_exp, h=6, w=8)

g4_MF_Pcla_exp <- goplot(go_enrich_MF_Pcla_exp, showCategory = 10)
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/EncrichedGOGraph_MF.pdf", g4_MF_Pcla_exp, h=6, w=8)

g6_MF_Pcla_exp <- cnetplot(go_enrich_MF_Pcla_exp, categorySize="pvalue")
ggplot2::ggsave("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/CategoryNet_MF.pdf", g6_MF_Pcla_exp, h=12, w=10)



##### plot #####

df_go_Pcla_exp <- list.files(paste0("../../data/analyzed_data/CAFE/results_cafe5_noPvirCopi/enrichment/clusterProfiler_Pcla_exp/"), 
                             pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~ data.table::fread(.)) %>%
  mutate(Description = str_to_sentence(Description),
         GeneRatio = str_split(GeneRatio, "/") %>% 
           map_chr(1) %>% 
           as.numeric() / str_split(GeneRatio, "/") %>% 
           map_chr(2) %>% 
           as.numeric()) %>%
  # group_by(Type) %>%
  arrange(qvalue) %>%
  slice_head(n = 10)

g4d <- ggplot(df_go_Pcla_exp,
              aes(y = reorder(Description, -qvalue), x = GeneRatio, col = qvalue)) +
  geom_point(shape = 16, size = 4) +#aes(size = Count),  , size= 4) +
  # scale_x_discrete(labels = scales::label_wrap(40)) +
  tidytext::scale_y_reordered() +
  scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
  # scale_shape_manual(values = 15, labels = "Molecular Function") +
  xlab("Gene ratio") +
  ylab("GO term") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())
g4d
ggplot2::ggsave("../../figures/tmp/Figure4d1.pdf", g4d, width = 4.5, height = 3)

