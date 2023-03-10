#### load libraries ####
library(tidyverse)
# install.packages("ggrepel")
library(ggrepel)
library(ape)
# install.packages("geiger")
library(geiger)
library(nlme)
# install.packages("phytools")
library(phytools)

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
df_prop <- read.delim("../../manuscript_tmp/data/propagule_size/results_20230220.txt") %>%
  mutate(Group = if_else(Group == 1, "Non-invasive", "Invasive") %>%
           factor(., levels = c("Invasive", "Non-invasive")),
         Decapods = if_else(Decapods == 1, "Decapods", "Others") %>%
           factor(., levels = c("Decapods", "Others")),
         spe_abb = paste0(substr(Species, 1, 1), ". ", str_split(Species, " ") %>% map_chr(2)))
colnames(df_prop) <- c("Invasive", "Decapods", "species", "phyname", "PD", "propsize", "spe_abb")

#### Figure 5a PD against propagule size ####
g5_1 <- ggplot(df_prop, aes(x = log10(propsize), y = PD)) +
  geom_smooth(method = 'lm', formula= y ~ x, color = "grey") +
  geom_point(aes(shape = Decapods, color = Invasive), size = 5) +
  geom_text_repel(aes(label = spe_abb),  fontface = "italic") +
  scale_color_manual(values = c("darkred", "grey")) +
  scale_shape_manual(values = c(17, 16)) +
  xlab(expression(Log[10]("Propagule size"))) +
  ylab(expression(italic(P)[D])) +
  theme_bw() +
  theme(legend.title = element_blank())

g5_1

# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/figures/tmp/Figure5.pdf", g5, w = 8, h = 8)


#### Figure 5b residuals ####
##### all invertebrates #####
#fit model
model <- lm(PD ~ log10(propsize), data=df_prop)
#view model summary
# summary(model)

df_prop2 <- df_prop %>%
  bind_cols(data.frame(res_std = rstandard(model)))

g5_2_all <- ggplot(df_prop2, aes(x = res_std, y = reorder(species, res_std))) +
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
allTree <- read.tree("../../manuscript_tmp/data/propagule_size/phylogenetic_tree/run.nex.all.treefile")
phy <- multi2di(allTree) #make rooted tree
plot(phy)
# # Calculate PICs
df_all_pgls <- df_prop %>%
  tibble::column_to_rownames(var = "phyname")
# pd_pic <- pic(df_all_pgls$PD, phy)
# prop_pic <- pic(df_all_pgls$propsize, phy)
# Make a model
# picModel <- lm(pd_pic ~ prop_pic - 1)
# Yes, significant
# summary(picModel)
# plot(pd_pic ~ prop_pic)
# abline(a = 0, b = coef(picModel))

pglsModel <- gls(PD ~ propsize, correlation = corBrownian(phy = allTree),
                 data = df_all_pgls, method = "ML")
summary(pglsModel)
plot(df_all_pgls$PD ~ df_all_pgls$propsize)
abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])

# ###### decapods ######
# model_dec <- lm(PD ~ log10(propsize), data = df_prop %>% filter(Decapods == "Decapods"))
# summary(model_dec)
# df_prop2_dec <- df_prop %>%
#   filter(Decapods == "Decapods") %>%
#   bind_cols(data.frame(res_std = rstandard(model_dec)))
# 
# g5_2_dec <- ggplot(df_prop2_dec, aes(x = res_std, y = reorder(species, res_std))) +
#   geom_bar(stat = "identity", aes(fill = Invasive)) +
#   # geom_text_repel(aes(label = species), fontface = "italic") + #nudge_x = 0.1, 
#   xlab(expression(paste("Residuals of ", {italic(P)[D]}, " against propagule size", sep = " "))) +
#   scale_fill_manual(values = c("darkred", "grey")) +
#   theme_classic() +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_text(face = "italic"),#element_blank(),
#         axis.title = element_blank(),
#         legend.position = "none")
# 
# g5_2 <- g5_2_all +
#   ggplot2::annotation_custom(grob = ggplotGrob(g5_2_dec), xmin = -2, xmax = 3, ymin = 1, ymax = 15)

# g5_2


g5 <- ggpubr::ggarrange(g5_1, g5_2_all, ncol = 2, align = "h")
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/figures/tmp/Figure5.pdf", g5, w = 14, h = 6)


# ###### PGLS analysis ######
# # Phylogenetic tree
# dacapodTree <- read.tree("../../ザリガニ論文_Sato_Makino/data/cafe/withoutPvir/raw_data/phylogenetic_tree/ultrametric_iqtree_name.nwk")
# plot(dacapodTree)
# # Calculate PICs
# df_decapod_pgls <- df_prop %>%
#   filter(Decapods == "Decapods") %>%
#   tibble::column_to_rownames(var = "species")
# pd_pic <- pic(df_decapod_pgls$PD, dacapodTree)
# prop_pic <- pic(df_decapod_pgls$propsize, dacapodTree)
# 
# # Make a model
# picModel <- lm(pd_pic ~ prop_pic - 1)
# 
# # Yes, significant
# summary(picModel)
# 
# plot(pd_pic ~ prop_pic)
# abline(a = 0, b = coef(picModel))




#### Figure 5c number of rapidly evolving genes ####
df_cafe_p <- read.table("../../ザリガニ論文_Sato_Makino/data/cafe/withoutPvir/raw_data/results_cafe5_noPvir/Base_branch_probabilities.txt", header = TRUE) %>%
  pivot_longer(cols = !GeneID, names_to = "node", values_to = "pval")
df_cafe_change <- bind_rows(read.table("../../ザリガニ論文_Sato_Makino/data/cafe/withoutPvir/raw_data/results_cafe5_noPvir/Base_change.txt", header = TRUE)) %>%
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
  filter(node == "Pcla.7", sign == "+")

#### Figure 5d enrichment analysis ####
list_cafe_all_genes <- read.table("../../ザリガニ論文_Sato_Makino/data/cafe/withoutPvir/raw_data/results_cafe5_noPvir/Base_change.txt", header = TRUE) %>%
  pull(GeneID)
list_cafe_Pcla_expanded_genes <- df_cafe_Pcla_exp %>%
  pull(GeneID)


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
write.table(df_go_enrich_BP_Pcla_exp, "../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/GO_BP_Pcla_exp.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_Pcla_exp <- upsetplot(go_enrich_BP_Pcla_exp)
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/UpsetPlot_BP.pdf", g1_BP_Pcla_exp, h=6, w=10)

g2_BP_Pcla_exp <- barplot(go_enrich_BP_Pcla_exp, 
                          drop = TRUE, 
                          showCategory = 10, 
                          title = "GO Biological Pathways",
                          font.size = 8)
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/BarPlot_BP.pdf", g2_BP_Pcla_exp, h=3, w=5)

g3_BP_Pcla_exp <- dotplot(go_enrich_BP_Pcla_exp, orderBy = "x")
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/DotPlot_BP.pdf", g3_BP_Pcla_exp, h=4, w=6)

g4_BP_Pcla_exp <- emapplot(go_enrich_BP_Pcla_exp %>% pairwise_termsim())
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/EncrichmentMap_BP.pdf", g4_BP_Pcla_exp, h=6, w=8)

g5_BP_Pcla_exp <- goplot(go_enrich_BP_Pcla_exp, showCategory = 10)
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/EncrichedGOGraph_BP.pdf", g5_BP_Pcla_exp, h=6, w=8)

g6_BP_Pcla_exp <- cnetplot(go_enrich_BP_Pcla_exp, categorySize="pvalue")
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/CategoryNet_BP.pdf", g6_BP_Pcla_exp, h=12, w=10)

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
# write.table(df_go_enrich_CC_Pcla_exp, "../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/GO_CC_Pcla_exp.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
# 
# g1_CC_Pcla_exp <- upsetplot(go_enrich_CC_Pcla_exp)
# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/UpsetPlot_CC.pdf", g1_CC_Pcla_exp, h=6, w=10)
# 
# g2_CC_Pcla_exp <- barplot(go_enrich_CC_Pcla_exp, 
#                           drop = TRUE, 
#                           showCategory = 10, 
#                           title = "GO Cellular Components",
#                           font.size = 8)
# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/BarPlot_CC.pdf", g2_CC_Pcla_exp, h=3, w=5)
# 
# g3_CC_Pcla_exp <- dotplot(go_enrich_CC_Pcla_exp, orderBy = "x")
# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/DotPlot_CC.pdf", g3_CC_Pcla_exp, h=4, w=5)
# 
# g4_CC_Pcla_exp <- emapplot(go_enrich_CC_Pcla_exp %>% pairwise_termsim())
# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/EncrichmentMap_CC.pdf", g4_CC_Pcla_exp, h=6, w=8)
# 
# g5_CC_Pcla_exp <- goplot(go_enrich_CC_Pcla_exp, showCategory = 10)
# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/EncrichedGOGraph_CC.pdf", g5_CC_Pcla_exp, h=6, w=8)
# 
# g6_CC_Pcla_exp <- cnetplot(go_enrich_CC_Pcla_exp, categorySize="pvalue")
# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/CategoryNet_CC.pdf", g6_CC_Pcla_exp, h=12, w=10)

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
write.table(df_go_enrich_MF_Pcla_exp, "../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/GO_MF_Pcla_exp.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_Pcla_exp <- upsetplot(go_enrich_MF_Pcla_exp)
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/UpsetPlot_MF.pdf", g1_MF_Pcla_exp, h=6, w=10)

g2_MF_Pcla_exp <- barplot(go_enrich_MF_Pcla_exp, 
                          drop = TRUE, 
                          showCategory = 10, 
                          title = "GO Molecular Functions",
                          font.size = 8)
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/BarPlot_MF.pdf", g2_MF_Pcla_exp, h=3, w=5)

g3_MF_Pcla_exp <- dotplot(go_enrich_MF_Pcla_exp, orderBy = "x")
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/DotPlot_MF.pdf", g3_MF_Pcla_exp, h=4, w=6)

g4_MF_Pcla_exp <- emapplot(go_enrich_MF_Pcla_exp %>% pairwise_termsim())
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/EncrichmentMap_MF.pdf", g4_MF_Pcla_exp, h=6, w=8)

g5_MF_Pcla_exp <- goplot(go_enrich_MF_Pcla_exp, showCategory = 10)
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/EncrichedGOGraph_MF.pdf", g5_MF_Pcla_exp, h=6, w=8)

g6_MF_Pcla_exp <- cnetplot(go_enrich_MF_Pcla_exp, categorySize="pvalue")
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/CategoryNet_MF.pdf", g6_MF_Pcla_exp, h=12, w=10)



##### plot #####

df_go_Pcla_exp <- list.files(paste0("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_Pcla_exp/"), 
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

g5d <- ggplot(df_go_Pcla_exp,
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
g5d
ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/figures/tmp/Figure5d1.pdf", g5d, width = 4.5, height = 3)

# g5d <- cnetplot(go_enrich_MF_Pcla_exp, 
#                 categorySize="qvalue",
#                 showCategory =8,
#                 vertex.label.font=6,
#                 layout = "kk",
#                 cex_category = 0.5,
#                 colorEdge= FALSE)
# 
# g5d  
# ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/figures/tmp/Figure5d2.pdf", g5d, width = 7, height = 6)

# write.table(list_cafe_all_genes, "../../ザリガニ論文_Sato_Makino/data/cafe/withoutPvir/list_genes/list_background.txt", append = F, quote = F, row.names = F, col.names = F) 
# write.table(list_cafe_Pcla_expanded_genes, "../../ザリガニ論文_Sato_Makino/data/cafe/withoutPvir/list_genes/list_Pcla_expanded.txt", append = F, quote = F, row.names = F, col.names = F)

