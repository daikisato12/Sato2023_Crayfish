#### load libraries ####
library(tidyverse)
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
# install.packages("flashClust")
library(flashClust)

# BiocManager::install("edgeR")
# library(edgeR)
# library(ggpubr)
# library(iterators)
# library(gtools) #mixedsort

# install.packages("tidytext")
library(tidytext)


#### Figure 3a PCA ####
## unique mapping ##
df_tpm_unique <- read.table("../../data/analyzed_data/RNAseq/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104.tsv", header = TRUE) 
sample <- colnames(df_tpm_unique)[2:16]
df_tpm_unique2 <- df_tpm_unique
pop <- factor(str_sub(colnames(df_tpm_unique)[2:16], 1, 3))
days <- factor(str_split(colnames(df_tpm_unique)[2:16], "_") %>% map_chr(1) %>% parse_number())
sampleinfo <- data.frame(name = colnames(df_tpm_unique)[2:16],
                         pop = pop,
                         days = days) %>%
  mutate(group = paste0(pop, "_", days))

pcDat <- prcomp(t(df_tpm_unique2 %>% 
                    dplyr::select(!c(ID))), 
                scale = TRUE)

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
list_degs_sen0_sen7 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05)
nrow(list_degs_sen0_sen7)

list_degs_sen0_sap0 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SapporoDay0_qval < 0.05)
nrow(list_degs_sen0_sap0)

list_degs_sen0_sap7 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SapporoDay7_qval < 0.05)
nrow(list_degs_sen0_sap7)

list_degs_sen0_sap31 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SapporoDay31_qval < 0.05)
nrow(list_degs_sen0_sap31)

list_degs_sen7_sap0 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay7_SapporoDay0_qval < 0.05)
nrow(list_degs_sen7_sap0)

list_degs_sen7_sap7 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay7_SapporoDay7_qval < 0.05)
nrow(list_degs_sen7_sap7)

list_degs_sen7_sap31 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay7_SapporoDay31_qval < 0.05)
nrow(list_degs_sen7_sap31)

list_degs_sap0_sap7 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SapporoDay0_SapporoDay7_qval < 0.05)
nrow(list_degs_sap0_sap7)

list_degs_sap0_sap31 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SapporoDay0_SapporoDay31_qval < 0.05)
nrow(list_degs_sap0_sap31)

list_degs_sap7_sap31 <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SapporoDay7_SapporoDay31_qval < 0.05)
nrow(list_degs_sap7_sap31)


#### Figure 3c WGCNA ####
df_tpm <- read.delim("../../data/analyzed_data/RNAseq/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104.tsv", header = TRUE, sep = "\t") %>%
  dplyr::rename(GeneID = ID) %>%
  left_join(read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
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

list_degs <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_withinpop.tsv", header = TRUE, sep = "\t") %>%
  mutate(rowname = paste0(GeneName, "_", FlyID, "_", GeneID))

datExpr_t <- df_tpm2[list_degs$rowname,] %>%
  na.omit()
# write.table(datExpr_t, "../figures/tmp/Table_WGCNA_most_variable_1000genes.tsv", sep="\t", row.names = T, quote = F, col.names = T)
write.table(datExpr_t, "../../data/analyzed_data/RNAseq/WGCNA/DEG_WGCNA.tsv", sep="\t", row.names = T, quote = F)

SubGeneNames <- rownames(datExpr_t)

datExpr <- t(log2(datExpr_t+1))

powers <- c(1:50)
sft <- pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers,corFnc = cor, corOptions = list(use = 'p'), networkType = "signed")

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
diag(dissTOM) <- NA

module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m <- t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))

G.expression <- datExpr
G.expressionColor <- numbers2colors(G.expression, signed = F, colors = gplots::colorpanel(100, low = "#F6F6C9", mid = "#BAD1C2", high = "#A3C7D6"))#gray.colors(100)) #colorpanel(100,low = "#839b5c", high = "#74325c", mid = "#E7E6D5")

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

write.table(tt2, "../../data/analyzed_data/RNAseq/WGCNA/WGCNA_Module.tsv", sep="\t", row.names = F, quote = F)


#### Figure 3d Enrichment analysis ####
##### write genes #####
# list of Sendai-specific DEGs
list_degs_sendai <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval),
         is.na(SapporoDay0_SapporoDay31_qval))
nrow(list_degs_sendai)
write.table(list_degs_sendai, "../../data/analyzed_data/RNAseq/genes/DEG_SendaiSpecific.tsv", sep="\t", row.names = F, quote = F)

# list of Sapporo-specific DEGs
list_degs_sapporo <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         SapporoDay0_SapporoDay7_qval < 0.05 | SapporoDay0_SapporoDay31_qval  < 0.05)
nrow(list_degs_sapporo)
write.table(list_degs_sapporo, "../../data/analyzed_data/RNAseq/genes/DEG_SapporoSpecific.tsv", sep="\t", row.names = F, quote = F)

# list of Population-shared DEGs
list_degs_shared <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_all.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval),
         SapporoDay0_SapporoDay31_qval  < 0.05)
nrow(list_degs_shared)
write.table(list_degs_shared, "../../data/analyzed_data/RNAseq/genes/DEG_Shared.tsv", sep="\t", row.names = F, quote = F)


##### run Figure3_enrichment.R #####

##### load results of enrichment analysis #####
df_go1 <- list.files(paste0("../../data/analyzed_data/RNAseq/enrichment/SapporoSpecificDEGs/"), 
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
  slice_head(n = 12) %>%
  mutate(type = "Sapporo-specific")

df_go2 <- list.files(paste0("../../data/analyzed_data/RNAseq/enrichment/SendaiSpecificDEGs/"), 
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
  slice_head(n = 12) %>%
  mutate(type = "Sendai-specific")

df_go3 <- list.files(paste0("../../data/analyzed_data/RNAseq/enrichment/SharedDEGs/"), 
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
  slice_head(n = 12) %>%
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

ggplot2::ggsave("../figures/tmp/Figure3d_GOenrichment.pdf", g3d, width = 6, height = 3.5)
