#### load libraries ####
library(tidyverse)
# install.packages("ggrepel")
library(ggrepel)
# install.packages("devtools", dependencies = TRUE)
# library(devtools)
# install_github('sinhrks/ggfortify')
# library(ggfortify)
library(data.table) #fread
library(scales) #label_wrap
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("WGCNA")
library(WGCNA)
library(flashClust)

# BiocManager::install("edgeR")
# library(edgeR)
# library(ggpubr)
# library(iterators)
# library(gtools) #mixedsort

# install.packages("tidytext")
library(tidytext)

setwd("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/codes/")

#### Figure 3a PCA ####
## unique mapping ##
df_tpm_unique <- read.table("../../manuscript_tmp/data/RNA-seq/20230104/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104.tsv", header = TRUE) 
# df_tpm_unique %>%
#   mutate(numTPM1 = rowSums(.[sample] > 1)) %>%
#   filter(numTPM1 > 1) %>%
#   dplyr::select(!numTPM1) %>%
#   write.table("../../manuscript_tmp/data/RNA-seq/20221230/DEG20221230/crayfish_hisat2_stringtie_uniq_tpm_20221230.tsv", row.names = F, quote = F, sep = "\t")
# df_tpm_unique %>%
#   mutate(numTPM1 = rowSums(.[sample] > 1)) %>%
#   filter(numTPM1 > 1, str_detect(ID, "^LOC")) %>%
#   dplyr::select(!numTPM1) %>%
#   write.table("../../manuscript_tmp/data/RNA-seq/20221230/DEG_annotated_only/crayfish_hisat2_stringtie_uniq_annotated_tpm_20221230.tsv", row.names = F, quote = F, sep = "\t")
sample <- colnames(df_tpm_unique)[2:16]
df_tpm_unique2 <- df_tpm_unique #%>%
  # mutate(numTPM1 = rowSums(.[sample] > 1)) %>%
  # filter(numTPM1 > 1) %>%
  # dplyr::select(!numTPM1)
pop <- factor(str_sub(colnames(df_tpm_unique)[2:16], 1, 3))
days <- factor(str_split(colnames(df_tpm_unique)[2:16], "_") %>% map_chr(1) %>% parse_number())
sampleinfo <- data.frame(name = colnames(df_tpm_unique)[2:16],
                         pop = pop,
                         days = days) %>%
  mutate(group = paste0(pop, "_", days))

# pop <- factor(str_split(colnames(df_tpm_unique2)[2:16], "_") %>% map_chr(1))
# days <- factor(str_split(colnames(df_tpm_unique2)[2:16], "_") %>% map_chr(2))
# sampleinfo <- data.frame(name = colnames(df_tpm_unique2)[2:16],
#                          pop = pop,
#                          days = days) %>%
#   mutate(group = paste0(pop, "_", days))

pcDat <- prcomp(t(df_tpm_unique2 %>% 
                    dplyr::select(!c(ID))), 
                scale = TRUE)
# ggbiplot(pcDat, obs.scale = 1, var.scale = 1, groups = sampleinfo$group, 
#          ellipse = TRUE, circle = FALSE, varname.size=0, var.axes = F)
df_pca <- data.frame(PC1 = scale(as.data.frame(pcDat$x)$PC1),
                     PC2 = scale(as.data.frame(pcDat$x)$PC2),
                     pop = sampleinfo$pop,
                     days = sampleinfo$days,
                     name = colnames(df_tpm_unique2)[2:16]) %>%
  mutate(pop = if_else(pop == "SEN", "Sendai", "Sapporo"),
         days = case_when(days == 0 ~ "Day 0",
                          days == 7 ~ "Day 7",
                          TRUE ~ "Day 31")) %>%
  transform(pop = factor(pop, levels = c("Sendai", "Sapporo")),
            days = factor(days, levels = c("Day 0", "Day 7", "Day 31")))
PC1_lab <- summary(pcDat)$importance[2,1] * 100
PC2_lab <- summary(pcDat)$importance[2,2] * 100


# plot PCA
g3a <- ggplot(df_pca, aes(x = PC1, y = PC2, col = pop, shape = days), size = 5) +
  geom_point(aes(x = PC1, y = PC2, col = pop, shape = days), size = 5) +
  geom_polygon(aes(x = PC1, y = PC2, fill = pop), alpha = 0.2) +#, position = position_jitter(width=0.3,height=0.3)) +
  # geom_text_repel(aes(label = name)) +
  scale_shape_manual(values=c(16, 17, 15), labels = c("Day 0", "Day 7", "Day 31")) +
  scale_color_manual(values = c("#a99e93", "#ebe1a9"), labels = c("Sendai", "Sapporo")) +
  scale_fill_manual(values = c("#a99e93", "#ebe1a9"), labels = c("Sendai", "Sapporo")) +
  xlab(paste0("PC1 (", round(PC1_lab, 2), "%)")) + 
  ylab(paste0("PC2 (", round(PC2_lab, 2), "%)")) +
  theme_bw() + 
  theme(legend.title = element_blank())

g3a
ggsave("../figures/tmp/Figure3a_unique_PCA.pdf", g3a, width = 4, height = 3)

df_tpm_unique_cv <- df_tpm_unique %>%
  pivot_longer(cols = contains("_"), names_sep = "_", names_to = c("pop", "rep"), values_to = "TPM") %>%
  group_by(ID, pop) %>%
  mutate(cv = sd(TPM) / mean(TPM))

#### Figure 3b Number of DEGs ####
# Manually created

# list of SendaiDEGsSapporo07NonDEGs
list_degs_SendaiDEGsSapporo07NonDEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval))
# (SapporoDay0_SapporoDay31_qval  < 0.05 | SapporoDay7_SapporoDay31_qval  < 0.05))

write.table(list_degs_SendaiDEGsSapporo07NonDEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiDEGsSapporo07NonDEGs.tsv", sep="\t", row.names = F, quote = F)

# list of SendaiOnlyDEGs
list_degs_SendaiOnlyDEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval), is.na(SapporoDay0_SapporoDay31_qval))
# (SapporoDay0_SapporoDay31_qval  < 0.05 | SapporoDay7_SapporoDay31_qval  < 0.05))

write.table(list_degs_SendaiOnlyDEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiOnlyDEGs.tsv", sep="\t", row.names = F, quote = F)


# list of SendaiDEGsSapporo07NonDEGs31DEGs
# list_degs_SendaiDEGsSapporo07NonDEGs31DEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
#   filter(SendaiDay0_SendaiDay7_qval < 0.05, 
#          is.na(SapporoDay0_SapporoDay7_qval),
#          SapporoDay0_SapporoDay31_qval < 0.05 | SapporoDay7_SapporoDay31_qval < 0.05)
# # (SapporoDay0_SapporoDay31_qval  < 0.05 | SapporoDay7_SapporoDay31_qval  < 0.05))
# 
# write.table(list_degs_SendaiDEGsSapporo07NonDEGs31DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiDEGsSapporo07NonDEGs31DEGs.tsv", sep="\t", row.names = F, quote = F)

# list of SendaiDEGsSapporo07NonDEGs31DEGs
list_degs_SendaiDEGsSapporo07NonDEGs31DEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval),
         SapporoDay0_SapporoDay31_qval < 0.05)
# (SapporoDay0_SapporoDay31_qval  < 0.05 | SapporoDay7_SapporoDay31_qval  < 0.05))

write.table(list_degs_SendaiDEGsSapporo07NonDEGs31DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiDEGsSapporo07NonDEGs31DEGs.tsv", sep="\t", row.names = F, quote = F)


# list of SendaiDEGsorSapporo07NonDEGs031DEGs
list_degs_SendaiDEGsorSapporo07NonDEGs031DEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SapporoDay0_SapporoDay7_qval),
         SendaiDay0_SendaiDay7_qval < 0.05 |
         SapporoDay0_SapporoDay31_qval < 0.05)
# (SapporoDay0_SapporoDay31_qval  < 0.05 | SapporoDay7_SapporoDay31_qval  < 0.05))

write.table(list_degs_SendaiDEGsorSapporo07NonDEGs031DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiDEGsorSapporo07NonDEGs031DEGs.tsv", sep="\t", row.names = F, quote = F)


# list of SendaiNonDEGsSapporo07DEGs
list_degs_SendaiNonDEGsSapporo07DEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         SapporoDay0_SapporoDay7_qval < 0.05)

write.table(list_degs_SendaiNonDEGsSapporo07DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiNonDEGsSapporo07DEGs.tsv", sep="\t", row.names = F, quote = F)

# list of SendaiNonDEGsSapporo031DEGs
list_degs_SendaiNonDEGsSapporo031DEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         is.na(SapporoDay0_SapporoDay7_qval), SapporoDay0_SapporoDay31_qval  < 0.05)

write.table(list_degs_SendaiNonDEGsSapporo031DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiNonDEGsSapporo031DEGs.tsv", sep="\t", row.names = F, quote = F)

# list of SendaiNonDEGsSapporo07031DEGs
list_degs_SendaiNonDEGsSapporo07031DEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         SapporoDay0_SapporoDay7_qval < 0.05 | SapporoDay0_SapporoDay31_qval < 0.05)

write.table(list_degs_SendaiNonDEGsSapporo07031DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiNonDEGsSapporo07031DEGs.tsv", sep="\t", row.names = F, quote = F)

# list of SendaiNonDEGsSapporo031731DEGs
list_degs_SendaiNonDEGsSapporo031731DEGs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         is.na(SapporoDay0_SapporoDay7_qval), SapporoDay0_SapporoDay31_qval  < 0.05, SapporoDay7_SapporoDay31_qval  < 0.05)

write.table(list_degs_SendaiNonDEGsSapporo031731DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiNonDEGsSapporo031731DEGs.tsv", sep="\t", row.names = F, quote = F)

# list of PopSpe07DEGs
list_degs_PopSpe07DEGs <- bind_rows(list_degs_SendaiDEGsSapporo07NonDEGs, list_degs_SendaiNonDEGsSapporo07DEGs)

write.table(list_degs_PopSpe07DEGs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_PopSpe07DEGs.tsv", sep="\t", row.names = F, quote = F)


# list of Sendai specific DEGs
list_degs_sendai <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval))

write.table(list_degs_sendai, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiSpecificDEGs.tsv", sep="\t", row.names = F, quote = F)



# list of Sapporo specific DEGs
list_degs_sapporo <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         SapporoDay0_SapporoDay7_qval < 0.05 | SapporoDay0_SapporoDay31_qval  < 0.05)

write.table(list_degs_sapporo, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SapporoSpecificDEGs.tsv", sep="\t", row.names = F, quote = F)

# list of Sapporo Day31 specific DEGs
list_degs_sapporo_Day31 <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         is.na(SapporoDay0_SapporoDay7_qval), SapporoDay0_SapporoDay31_qval  < 0.05 | SapporoDay7_SapporoDay31_qval < 0.05)

write.table(list_degs_sapporo_Day31, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SapporoDay31SpecificDEGs.tsv", sep="\t", row.names = F, quote = F)


# list of Population specific DEGs
list_degs_psdegs <- bind_rows(list_degs_sendai, list_degs_sapporo)
write.table(list_degs_psdegs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_PopulationSpecificDEGs.tsv", sep="\t", row.names = F, quote = F)

# list of Shared DEGs
list_degs_shareddegs <- read.delim("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         SapporoDay0_SapporoDay7_qval < 0.05 | SapporoDay0_SapporoDay31_qval < 0.05)
write.table(list_degs_shareddegs, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SharedDEGs.tsv", sep="\t", row.names = F, quote = F)


#### Figure 3c WGCNA ####
df_tpm <- read.delim("../data/RNA-seq/unique_mapping/rawdata/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104.tsv", header = TRUE, sep = "\t") %>%
  dplyr::rename(GeneID = ID) %>%
  left_join(read.delim("../figures/tmp/TableS2.tsv", header = TRUE, sep = "\t") %>%
              dplyr::select(GeneID, FlyID, GeneName)) %>%
  pivot_longer(cols = contains("_"), values_to = "TPM") %>%
  separate(col = name, into = c("group", "no"), sep = "_") %>%
  mutate(pop = str_sub(group, 1, 3),
         group = case_when(group == "SEN0" ~ "Sendai Day 0",
                           group == "SEN7" ~ "Sendai Day 7",
                           group == "TON0" ~ "Sapporo Day 0",
                           group == "TON7" ~ "Sapporo Day 7",
                           group == "TON30" ~ "Sapporo Day 31"))

df_tpm2 <- df_tpm %>%
  mutate(sample = paste0(group, " rep", no),
         rowname = paste0(GeneName, "_", FlyID, "_", GeneID)) %>%
  filter(rowname != "NA_NA_NA") %>%
  dplyr::select(rowname, sample, TPM) %>%
  pivot_wider(names_from = sample, values_from = TPM) %>%
  arrange(rowname) %>%
  tibble::column_to_rownames()

# select 1000 most variable genes by cv (coefficient of variation)
# N <- 1000
# df_cv <- df_tpm %>% 
#   #  rowwise() %>%
#   mutate(mean = apply(., 1, mean),
#          sd = apply(., 1, sd),
#          cv = sd / mean)
# 
# df_cv2 <- df_cv %>% arrange(desc(cv))
# df_1000genes <- rownames(df_cv2)[1:N] %>% as.data.frame()
# colnames(df_1000genes) <- "Geneid"

# use DEGs instead
list_degs <- read.delim("../figures/tmp/TableS2_withinDEGs.tsv", header = TRUE, sep = "\t") %>%
  mutate(rowname = paste0(GeneName, "_", FlyID, "_", GeneID))

datExpr_t <- df_tpm2[list_degs$rowname,] %>%
  na.omit()
# write.table(datExpr_t, "../figures/tmp/Table_WGCNA_most_variable_1000genes.tsv", sep="\t", row.names = T, quote = F, col.names = T)
write.table(datExpr_t, "../figures/tmp/Table_WGCNA_DEGs.tsv", sep="\t", row.names = T, quote = F)

SubGeneNames <- rownames(datExpr_t)

datExpr <- t(log2(datExpr_t+1))

powers <- c(1:50)
sft <- pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers,corFnc = cor, corOptions = list(use = 'p'), networkType = "signed")
#sft=pickSoftThreshold(datExpr, powerVector=powers, networkType = "unsigned", RsquaredCut=0.8, verbose=5)
#sft$powerEstimate

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h = 0.80, col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower <- 30 #best threshold? (by the way, here we can not specify more than 30)

#calclute the adjacency matrix
adj <- adjacency(datExpr, type = "signed", power = softPower)

#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM <- TOMsimilarityFromExpr(datExpr, networkType = "signed", TOMType = "signed", power = softPower)

colnames(TOM) = rownames(TOM) <- SubGeneNames
dissTOM <- 1-TOM

#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree <- flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize <- 20;

# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#set the diagonal of the dissimilarity to NA 
diag(dissTOM) <- NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
# sizeGrWindow(7,7)
# TOMplot(dissTOM^4, geneTree, as.character(dynamicColors))
# 
# module_colors= setdiff(unique(dynamicColors), "grey")
# for (color in module_colors){
#   module=SubGeneNames[which(dynamicColors==color)]
#   write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
# }

module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m <- t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
# m2 <- m[sample_amy,]
# heatmap(t(m2),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])

G.expression <- datExpr
G.expressionColor <- numbers2colors(G.expression, signed = F, colors = gplots::colorpanel(100, low = "#F6F6C9", mid = "#BAD1C2", high = "#A3C7D6"))#gray.colors(100)) #colorpanel(100,low = "#839b5c", high = "#74325c", mid = "#E7E6D5")
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

dynamicColors2 <- data.frame(dynamicColors, G.expressionColor[1,], G.expressionColor[2,],
                             G.expressionColor[3,], G.expressionColor[4,], G.expressionColor[5,],
                             G.expressionColor[6,], G.expressionColor[7,], G.expressionColor[8,],
                             G.expressionColor[9,], G.expressionColor[10,], G.expressionColor[11,],
                             G.expressionColor[12,], G.expressionColor[13,], G.expressionColor[14,],
                             G.expressionColor[15,])
col_order <- c("Sendai Day0 rep1", "Sendai Day0 rep2", "Sendai Day0 rep3",
               "Sendai Day7 rep1", "Sendai Day7 rep2", "Sendai Day7 rep3",
               "Sapporo Day0 rep1", "Sapporo Day0 rep2", "Sapporo Day0 rep3", 
               "Sapporo Day7 rep1", "Sapporo Day7 rep2", "Sapporo Day7 rep3", 
               "Sapporo Day31 rep1", "Sapporo Day31 rep2", "Sapporo Day31 rep3")
pdf("../figures/tmp/WGCNA_tree_all.pdf", w=7, h=5)
plotDendroAndColors(geneTree, dynamicColors2, groupLabels = c("Dynamic Tree Cut", col_order), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()


# show color bar
library(heatmap3)
library(gplots)
pdf("../figures/tmp/WGCNA_colorlegend.pdf", w=7, h=5)
colByValue(G.expression, col = gplots::colorpanel(100, low = "#F6F6C9", mid = "#BAD1C2", high = "#A3C7D6"))
dev.off()


unique(dynamicColors2[geneTree$order,1])
# Module cyan/salmon/lightcyan/blue/green/black/greenyellow/brown/magenta/turquoise
tt <- datExpr_t[rownames(dynamicColors2[str_detect(dynamicColors, "cyan|salmon|lightcyan|blue|green|black|greenyellow|brown|magenta|turquoise"),]) %>% as.integer(),] %>%
  as.matrix()

tt2 <- log2(tt+1) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  separate(col = rowname, into = c("Genename", "FlyID", "GeneID"), sep = "_") %>%
  bind_cols(data.frame(Module = dynamicColors[str_detect(dynamicColors, "cyan|salmon|lightcyan|blue|green|black|greenyellow|brown|magenta|turquoise")])) %>%
  relocate(GeneID, FlyID, Genename, Module) %>%
  arrange(Genename, GeneID)

write.table(tt2, "../figures/tmp/WGCNA_Module.tsv", sep="\t", row.names = F, quote = F)


# tt_diff <- datExpr_t[rownames(dynamicColors2[str_detect(dynamicColors, "tan|red|purple"),]) %>% as.integer(),] %>%
#   as.matrix()
# tt2_diff <- log2(tt_diff+1) %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column() %>%
#   separate(col = rowname, into = c("Genename", "FlyID", "GeneID"), sep = "_") %>%
#   bind_cols(data.frame(Module = dynamicColors[str_detect(dynamicColors, "tan|red|purple")])) %>%
#   relocate(GeneID, FlyID, Genename, Module) %>%
#   arrange(Genename, GeneID)
# 
# write.table(tt2_diff, "../figures/tmp/WGCNA_Module_diff.tsv", sep="\t", row.names = F, quote = F)


datExpr_t3 <- datExpr_t %>%
  tibble::rownames_to_column() %>%
  separate(col = rowname, into = c("Name", "FlyID", "GeneID"), sep = "_") %>%
  bind_cols(data.frame(color = dynamicColors))
datExpr_t3_deg <- inner_join(datExpr_t3, list_degs_SendaiDEGsSapporo07NonDEGs)

df_deg_enrich <- data.frame()
for (color_i in unique(dynamicColors)){
  # color_i <- unique(dynamicColors)[1]
  datExpr_temp <- datExpr_t3 %>%
    filter(color == color_i)
  
  datExpr_temp_deg <- inner_join(datExpr_temp, list_degs_SendaiDEGsSapporo07NonDEGs)
  p <- fisher.test(matrix(c(nrow(datExpr_t3), nrow(datExpr_t3_deg), nrow(datExpr_temp), nrow(datExpr_temp_deg)), nrow = 2))$p
  odds.ratio <- fisher.test(matrix(c(nrow(datExpr_t3), nrow(datExpr_t3_deg), nrow(datExpr_temp), nrow(datExpr_temp_deg)), nrow = 2))$estimate
  df_deg_enrich <- bind_rows(df_deg_enrich,
                             data.frame(color = color_i,
                                        num_all = nrow(datExpr_temp),
                                        num_deg = nrow(datExpr_temp_deg),
                                        p.value = p,
                                        odds.ratio = odds.ratio
                                        ))
}
df_deg_enrich_sig <- df_deg_enrich %>%
  filter(p.value < 0.05)

# list of SendaiDEGsSapporo07NonDEGs enriched in module
datExpr_t3_deg_enrichmodule <- datExpr_t3_deg %>%
  filter(color %in% unique(df_deg_enrich_sig$color))

write.table(datExpr_t3_deg_enrichmodule, "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/TableS5_SendaiDEGsSapporo07NonDEGsEnrichedModule.tsv", sep="\t", row.names = F, quote = F)


datExpr_temp_deg <- inner_join(tt2, list_degs_SendaiDEGsSapporo07NonDEGs)

p <- fisher.test(matrix(c(nrow(datExpr_t3), nrow(datExpr_t3_deg), nrow(tt2), nrow(datExpr_temp_deg)), nrow = 2))$p
odds.ratio <- fisher.test(matrix(c(nrow(datExpr_t3), nrow(datExpr_t3_deg), nrow(datExpr_temp), nrow(datExpr_temp_deg)), nrow = 2))$estimate

list_degs_SendaiDEGsSapporo07NonDEGslist_degs_SendaiDEGsSapporo07NonDEGs %>%
  # filter(SapporoDay0_SapporoDay31_qval < 0.05 | SapporoDay7_SapporoDay31_qval < 0.05) %>%
  nrow()

# #### Figure 3c Enrichment analysis ####
# ### run Figure3_enrichment.R
# 
# df_go <- list.files(paste0("../data/enrichment/clusterProfiler_DEG2/"), 
#                     pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
#   map_df(~ data.table::fread(.))
# 
# g3c <- ggplot(df_go %>%
#                 mutate(Description = str_to_sentence(Description),
#                        GeneRatio = str_split(GeneRatio, "/") %>% map_chr(1) %>% as.numeric() / str_split(GeneRatio, "/") %>% map_chr(2) %>% as.numeric()) %>%
#                 # group_by(Type) %>%
#                 arrange(qvalue) %>%
#                 slice_head(n = 20), 
#               aes(x = GeneRatio, y = reorder(Description, GeneRatio), col = qvalue, size = Count)) +
#   geom_point(aes(shape = Type)) +
#   scale_y_discrete(labels = scales::label_wrap(40)) +
#   scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
#   scale_shape_manual(values = c(16, 17, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
#   xlab("Gene ratio") +
#   ylab("GO term") +
#   # facet_wrap(~ Type, nrow = 3, scales = "free_y") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# g3c
# ggsave("../figures/tmp/Figure3c_GOenrichment.pdf", g3c, width = 6, height = 5)

#### Figure 3d Enrichment analysis ####
### run Figure3_enrichment.R

df_go1 <- list.files(paste0("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_SapporoDEGs/"), 
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
  slice_head(n = 14) %>%
  mutate(type = "Sapporo-specific")

# g3d1 <- ggplot(df_go1, aes(x = GeneRatio, y = reorder(Description, -qvalue), col = qvalue)) +
#   geom_point(aes(shape = Type), size= 4) +
#   scale_y_discrete(labels = scales::label_wrap(40)) +
#   scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
#   scale_shape_manual(values = c(16, 17, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
#   xlab("Gene ratio") +
#   ylab("GO term") +
#   # facet_wrap(~ Type, nrow = 3, scales = "free_y") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# g3d1

df_go2 <- list.files(paste0("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_SendaiOnlyDEGs/"), 
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
  slice_head(n = 10) %>%
  mutate(type = "Sendai-specific")

# g3d2 <- ggplot(df_go2, aes(x = GeneRatio, y = reorder(Description, -qvalue), col = qvalue)) +
#   geom_point(aes(shape = Type), size= 4) +
#   scale_y_discrete(labels = scales::label_wrap(40)) +
#   scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
#   scale_shape_manual(values = c(16, 17, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
#   xlab("Gene ratio") +
#   ylab("GO term") +
#   # facet_wrap(~ Type, nrow = 3, scales = "free_y") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# g3d2
# 
# g3d <- ggpubr::ggarrange(g3d1, g3d2, nrow = 2, align = "v")


df_go3 <- list.files(paste0("../../ザリガニ論文_Sato_Makino/data/enrichment/clusterProfiler_SendaiDEGsSapporo07NonDEGs31DEGs/"), 
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
  slice_head(n = 15) %>%
  mutate(type = "Shared")

g3d <- ggplot(bind_rows(df_go1, df_go2) %>%
                bind_rows(df_go3) %>%
                filter(Type %in% c("BP", "MF")), 
              aes(y = reorder_within(Description, -qvalue, type), x = GeneRatio, col = qvalue)) +
  geom_point(aes(shape = Type, size = Count)) +#, size= 4) +
  # scale_x_discrete(labels = scales::label_wrap(40)) +
  scale_y_reordered() +
  scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
  # scale_shape_manual(values = c(16, 17, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
  scale_shape_manual(values = c(17, 16), labels = c("Biological process", "Molecular Function")) +
  xlab("Gene ratio") +
  ylab("GO term") +
  facet_grid(type ~ ., scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())

g3d

ggplot2::ggsave("../../ザリガニ論文_Sato_Makino/figures/tmp/Figure3d_GOenrichment.pdf", g3d, width = 6, height = 3.5)
