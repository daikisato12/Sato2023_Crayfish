#### WGCNA ####
library(tidyverse)
# BiocManager::install("WGCNA")
library(WGCNA)
library(flashClust)
# install.packages("gplots")
library(gplots)

setwd("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/codes/")

df_tpm <- read.table("../data/RNA-seq/rawdata/gene_level_tpm/crayfish_hisat2_stringtie_gene_tpm20221213_LOC.tsv", header = TRUE) %>%
  left_join(read.table("../figures/tmp/TableS2.tsv", header = TRUE) %>%
              dplyr::select(GeneID, FlyID, GeneName)) %>%
  pivot_longer(cols = contains("_"), values_to = "TPM") %>%
  separate(col = name, into = c("group", "no"), sep = "_") %>%
  mutate(pop = str_sub(group, 1, 3),
         group = case_when(group == "SEN0" ~ "Sendai Day 0",
                           group == "SEN7" ~ "Sendai Day 7",
                           group == "TON0" ~ "Sapporo Day 0",
                           group == "TON7" ~ "Sapporo Day 7",
                           group == "TON30" ~ "Sapporo Day 30"))

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
list_degs <- read.table("../figures/tmp/TableS2_withinDEGs.tsv", header = TRUE) %>%
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
pdf("../figures/tmp/WGCNA_tree_all.pdf", w=7, h=5)
plotDendroAndColors(geneTree, dynamicColors2, groupLabels = c("Dynamic Tree Cut", col_order), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# show color bar
library(heatmap3)
library(gplots)
pdf("../figures/tmp/WGCNA_colorlegend.pdf", w=7, h=5)
colByValue(G.expression, col = gplots::colorpanel(100, low = "#F6F6C9", mid = "#BAD1C2", high = "#A3C7D6"))
dev.off()


#### Module genes ####
#### Module black ####
tt <- datExpr_t[rownames(dynamicColors2[dynamicColors=="black",]) %>% as.integer(),] %>%
  as.matrix()

tt2 <- log2(tt+1) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  separate(col = rowname, into = c("Genename", "FlyID", "GeneID"), sep = "_")

write.table(tt2, "../figures/tmp/WGCNA_Module_black.tsv", sep="\t", row.names = F, quote = F)

#### Module red ####
tt <- datExpr_t[rownames(dynamicColors2[dynamicColors=="red",]) %>% as.integer(),] %>%
  as.matrix()

tt2 <- log2(tt+1) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  separate(col = rowname, into = c("Genename", "FlyID", "GeneID"), sep = "_")

write.table(tt2, "../figures/tmp/WGCNA_Module_red.tsv", sep="\t", row.names = F, quote = F)

#### Module pink ####
tt <- datExpr_t[rownames(dynamicColors2[dynamicColors=="pink",]) %>% as.integer(),] %>%
  as.matrix()

tt2 <- log2(tt+1) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  separate(col = rowname, into = c("Genename", "FlyID", "GeneID"), sep = "_")

write.table(tt2, "../figures/tmp/WGCNA_Module_pink.tsv", sep="\t", row.names = F, quote = F)


#### Module turquoise/black/red/pink/brown/purple/yellow ####
tt <- datExpr_t[rownames(dynamicColors2[str_detect(dynamicColors, "turquoise|black|red|pink|brown|purple|yellow"),]) %>% as.integer(),] %>%
  as.matrix()

tt2 <- log2(tt+1) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  separate(col = rowname, into = c("Genename", "FlyID", "GeneID"), sep = "_")

write.table(tt2, "../figures/tmp/WGCNA_Module.tsv", sep="\t", row.names = F, quote = F)
