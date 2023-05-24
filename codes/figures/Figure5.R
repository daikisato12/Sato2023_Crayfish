#### load libraries ####
library(tidyverse)
library(patchwork)

#### Figure 5a Manhattan plot ####
df_pbs <- read.table("../../data/analyzed_data/PBS/HAF_PBS_5kbp_step1kbp.tsv", header = TRUE) %>%
  pivot_longer(cols = !c(Chr, Coord, NSNPs), names_sep = "\\.", names_to = c("TargetPop", "comb1", "comb2", "comb3"), values_to = "PBS") %>%
  unite(col = Combination, starts_with("comb"), sep = "-") %>%
  group_by(Combination, TargetPop) %>%
  mutate(threshold1 = quantile(na.omit(PBS),0.99),
         threshold0.1 = quantile(na.omit(PBS),0.999))

unique(df_pbs$Chr)
li <- seq(1, length(unique(df_pbs$Chr)))
names(li) <- unique(df_pbs$Chr)

df_pbs2 <- df_pbs %>%
  mutate(CHR_tmp = li[Chr]) %>%
  dplyr::rename(BP = Coord)

don_PBS <- df_pbs2 %>%
  
  # filter(df_PBS != "na") %>%
  
  # Compute CHR_tmpomosome size
  group_by(CHR_tmp) %>% 
  summarise(CHR_tmp_len = max(BP)) %>% 
  
  # Calculate cumulative position of each CHR_tmpomosome
  mutate(tot = cumsum(as.numeric(CHR_tmp_len))-CHR_tmp_len) %>%
  dplyr::select(-CHR_tmp_len) %>%
  
  # Add this info to the initial dataset
  left_join(df_pbs2, ., by=c("CHR_tmp"="CHR_tmp")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR_tmp, BP) %>%
  mutate(BPcum=BP+tot) %>%
  mutate(sig_color = case_when(as.numeric(CHR_tmp) %% 2 == 1 ~ "black",
                               as.numeric(CHR_tmp) %% 2 == 0 ~ "grey"))

axisdf <- don_PBS %>% 
  group_by(CHR_tmp) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g5a <- ggplot(don_PBS %>%
                dplyr::rename(POP = Combination) %>%
                filter(TargetPop == "Sapporo"), 
              aes(x=BPcum, y=as.numeric(PBS), color = sig_color)) + # color= as.factor(CHR_tmp) # alpha = as.factor(CHR_tmp %% 2), color = POP)) +
  
  # Show all points
  geom_point(size = 1.3, alpha = 0.4, shape = 16) +
  # scale_alpha_manual(values = c(0.5, 1)) +
  # scale_color_viridis_d() +
  scale_color_manual(values = c("black" = "#d9d9b4",
                                "grey" = "#fff1cf")) +
  # scale_color_manual(values = rep(c("grey", "black"), max(don_PBS$CHR_tmp))) +
  #  scale_size_manual(values = c(1.3, 2.4)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  #  scale_fill_viridis_d() +
  # scale_shape_manual(values = c(16, 17, 15)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  # geom_hline(aes(yintercept = threshold1), linewidth = 1, col = "black", linetype = "dotted") +
  geom_hline(aes(yintercept = threshold0.1), linewidth = 1, col = "darkred", linetype = "dotted") +
  # geom_hline(yintercept = quantile(na.omit(don_PBS$PBS),0.99), linewidth = 1, col = "black", linetype = "dotted") +
  # geom_hline(yintercept = quantile(na.omit(don_PBS$PBS),0.999), linewidth = 1, col = "darkred", linetype = "dotted") +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR_tmp, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  #  scale_y_continuous(breaks=seq(0,1,by=0.2)) +     # remove space between plot area and x axis
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Linkage Group") +
  ylab("PBS (Sapporo)") +
  coord_cartesian(ylim = c(0,NA)) +
  # Custom the theme:
  facet_wrap(~ POP, nrow = 3) +
  theme_bw() +
  theme( 
    # axis.text.x = element_blank(),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank(),
  ) 
# g5a
ggsave("../figures/tmp/Figure5a_Manhattan_plot_5kbp_step1kbp.png", g5a, w=6, h=3)

don_PBS_5kbp_1kbp_0.1 <- don_PBS %>%
  filter(PBS > threshold0.1, TargetPop == "Sapporo") %>%
  dplyr::select(!c(CHR_tmp, tot, BPcum, sig_color)) %>%
  pivot_wider(id_cols = c(Chr, BP, NSNPs), names_from = Combination, values_from = PBS) %>%
  write.table("../../data/analyzed_data/PBS/HAF_PBS_5kbp_step1kbp_top0.1.tsv",row.names = F, quote = F, sep = "\t")

comb <- c("Kamakura-Okinawa-Sapporo",
          "Aomori-Kamakura-Sapporo",
          "Aomori-Okinawa-Sapporo")

don_PBS_5kbp_1kbp_0.1_2comb <- don_PBS %>%
  filter(PBS > threshold0.1, TargetPop == "Sapporo") %>%
  dplyr::select(!c(CHR_tmp, tot, BPcum, sig_color)) %>%
  pivot_wider(id_cols = c(Chr, BP, NSNPs), names_from = Combination, values_from = PBS) %>%
  mutate(numTop = rowSums(!is.na(.[comb]))) %>%
  filter(numTop > 1) %>%
  dplyr::select(!numTop) %>%
  write.table("../../data/analyzed_data/PBS/HAF_PBS_5kbp_step1kbp-top0.1_2combsharedloci.tsv",row.names = F, quote = F, sep = "\t")

#### Figure 5b Number of genes ####
list_degs_DEGs <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_withinpop.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05 |
           SapporoDay0_SapporoDay7_qval < 0.05 | 
           SapporoDay0_SapporoDay31_qval  < 0.05 |
           SapporoDay7_SapporoDay31_qval  < 0.05) %>%
  pull(GeneID)

list_degs_sendai <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_withinpop.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval), 
         is.na(SapporoDay0_SapporoDay31_qval)) %>%
  pull(GeneID)

list_degs_sapporo <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_withinpop.tsv", header = TRUE, sep = "\t") %>%
  filter(is.na(SendaiDay0_SendaiDay7_qval), 
         SapporoDay0_SapporoDay7_qval < 0.05 | SapporoDay0_SapporoDay31_qval  < 0.05) %>%
  pull(GeneID)

list_degs_shared <- read.delim("../../data/analyzed_data/RNAseq/genes/DEG_withinpop.tsv", header = TRUE, sep = "\t") %>%
  filter(SendaiDay0_SendaiDay7_qval < 0.05, 
         is.na(SapporoDay0_SapporoDay7_qval),
         SapporoDay0_SapporoDay31_qval < 0.05) %>%
pull(GeneID)

list_selection <- read.delim("../../data/analyzed_data/PBS/HAF_PBS_5kbp_step1kbp-top0.1_2combsharedloci_gene.tsv") %>%
  filter(Chr == "LG03_CM035475.1") %>%
  separate_rows(Genes, sep = ", ") %>%
  pull(Genes) %>%
  na.omit() %>%
  unique()

list_selection <- read.delim("../../data/analyzed_data/PBS/HAF_PBS_5kbp_step1kbp-top0.1_2combsharedloci_gene.tsv") %>%
  # filter(Chr == "LG03_CM035475.1") %>%
  separate_rows(Genes, sep = ", ") %>%
  pull(Genes) %>%
  na.omit() %>%
  unique()

intersect(list_degs_sendai, list_selection)
intersect(list_degs_sapporo, list_selection)
intersect(list_degs_shared, list_selection)
intersect(list_degs_DEGs, list_selection)

setdiff(list_degs_DEGs, 
        c(list_degs_sendai,
          list_degs_sapporo,
          list_degs_shared))

#### Figure 5c Heatmap of allele frequencies ####
label_plus <- function(x){
  ifelse(sign(x) == -1, abs(x), x)
}

df_gtf <- read_tsv("../../data/analyzed_data/Ref/GCF_020424385.1_ASM2042438v2_genomic.gtf", comment = "#", col_names = F) %>%
  separate(X9, sep = ";", into = c("geneID", "transcriptID", "db_xref", "experiment", "gene", "model_evidence", "product", "transcript_biotype", "exon_number")) %>%
  dplyr::select(1:10) %>%
  mutate(geneID = str_extract(geneID, '\"[^()]+\"') %>%
           str_remove_all('\"'),
         transcriptID = str_extract(transcriptID, '\"[^()]+\"') %>%
           str_remove_all('\"'))
genename_list <- df_gtf[!is.na(df_gtf$transcriptID),]$geneID
names(genename_list) <- df_gtf[!is.na(df_gtf$transcriptID),]$transcriptID

df_pbs <- read.table("../../data/analyzed_data/PBS/HAF_PBS_5kbp_step1kbp.tsv", header = TRUE) %>%
  pivot_longer(cols = !c(Chr, Coord, NSNPs), names_sep = "\\.", names_to = c("TargetPop", "Comb1", "Comb2", "Comb3"), values_to = "PBS") %>%
  unite(col = Combination, starts_with("comb"), sep = "-")

df_g5c1 <- read.table("../../data/analyzed_data/PBS/merge_result.tsv", sep = "\t", header = TRUE) %>% 
  mutate(SNPID = row_number(),
         check0 = rowSums(.[3:9] == 0),
         check1 = rowSums(.[3:9] == 1)) %>%
  filter(check0 != 7 & check1 != 7) %>%
  dplyr::rename(Atchafalaya = atchafalaya2,
                Kamakura = kamakura2,
                Okinawa = okinawa2,
                Tonden = tonden2,
                Triunfo = triunfo2) %>%
  pivot_longer(cols = c("Aomori", "Atchafalaya", "Kamakura", "NewOrleans", "Okinawa", "Tonden", "Triunfo"), values_to = "AF") %>%
  mutate(AF = as.numeric(AF),
         name = case_when(name == "Tonden" ~ "Sapporo",
                          name == "NewOrleans" ~ "New Orleans",
                          TRUE ~ name)) %>%
  transform(name = factor(name, levels = c("Okinawa", "Aomori", "Sapporo", "Kamakura", "Atchafalaya", "Triunfo", "New Orleans")))

region <- c(22000000, 27000000)

df_g5c1_t <- df_g5c1 %>%
  filter(Chr == "LG03_CM035475.1", Pos > region[1], Pos < region[2]) %>%
  dplyr::rename(Coord = Pos)

df_g5c1_2 <- df_g5c1_t %>% #df_g5c1 %>%
  pivot_wider(id_cols = Coord, names_from = name, values_from = AF) %>%
  tibble::column_to_rownames(var = "Coord") %>%
  t()
df_g5c1_2 <- df_g5c1_2[c("Okinawa", "Aomori", "Sapporo", "Kamakura", "Atchafalaya", "Triunfo", "New Orleans"),]
h <- hclust(dist(df_g5c1_2))
d <- as.dendrogram(h)
dd <- ggdendro::dendro_data(d, type = "rectangle")

df_gtf_exon <- df_gtf %>%
  filter(X4 > region[1], X5 < region[2], X3 == "exon") %>%
  filter(X1 == "NC_059573.1")
xmin <- df_gtf_exon$X4
xmax <- df_gtf_exon$X5
df_gtf_gene <- df_gtf %>%
  filter(X4 > region[1], X5 < region[2], X3 == "gene") %>%
  filter(X1 == "NC_059573.1")

g5c1_dend <- ggplot(ggdendro::segment(dd)) +
  geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
  scale_x_continuous(labels = label_plus) +
  xlab("Height") +
  # geom_text(data = dd$labels, aes(x = -y, y = x, label = label), hjust = 0) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank())
g5c1_dend

g5c1_pbs <- ggplot(df_pbs %>%
                     dplyr::rename(CHR = Chr,
                                   BP = Coord) %>%
                     filter(CHR == "LG03_CM035475.1", 
                            BP > region[1], BP < region[2],
                            TargetPop == "Sapporo")) +
  annotate("rect", xmin = 24500000, xmax = 25500000, ymin = -Inf, ymax = Inf,  fill = "red", alpha=.2) +
  geom_point(aes(x = BP, y = PBS, color = Combination), alpha = 0.5, shape = 16, size = 0.5) +
  scale_color_viridis_d(labels = c("Aomori-Kamakura-Sapporo",
                                   "Aomori-Okinawa-Sapporo",
                                   "Kamakura-Okinawa-Sapporo")) +
  scale_x_continuous(labels = scales::label_comma(scale = 1/1000000)) + # 単位M(1000000)をつけるため値自体は1/1000000する
  coord_cartesian(xlim = c(23000000, 26000000), ylim = c(0, 1.8)) +
  xlab("Linkage group 3 (Mb)") +
  # geom_text(data = dd$labels, aes(x = -y, y = x, label = label), hjust = 0) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank())

g5c1_gene <- ggplot(df_gtf_gene %>%
                      filter(X4 > 24500000, X5 < 25500000), 
                    aes(x = X4, y = 1, label = geneID)) +
  annotate("rect", xmin = xmin, xmax = xmax, ymin = 0.8, ymax = 1.2, color = "black", fill = "black") +
  annotate("segment", 
           x = df_gtf_gene %>% pull(X5), 
           xend = df_gtf_gene %>% pull(X4), 
           y = 1, yend = 1) +
  annotate("segment", 
           x = df_gtf_gene %>% pull(X5), 
           xend = df_gtf_gene %>% pull(X4), 
           y = 1, yend = 1) + 
  # ggrepel::geom_label_repel(max.overlaps = 10) +
  annotate("text",
           size = 1,
           label = df_gtf_gene %>% pull(geneID),
           x = (df_gtf_gene %>% pull(X4) + df_gtf_gene %>% pull(X5) ) /2,
           y = c(rep(seq(0.5, 0.7, by = 0.1), as.integer(length(df_gtf_gene$X4)/3)), 0.5, 0.6)) +
  xlab("Linkage group 3 (Mb)") +
  # coord_cartesian(xlim = c(3124000, 3160000)) +
  coord_cartesian(xlim = c(24500000, 25500000)) +
  scale_x_continuous(labels = scales::label_comma(scale = 1/1000000)) + # 単位M(1000000)をつけるため値自体は1/1000000する
  # ggdendro::theme_dendro() +
  theme_classic() +
  theme(plot.margin = unit(c(0,0,0,0), "null"), #top, right, bottom, left 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_blank())
# g5c1_gene

g5c1 <- ggplot(df_g5c1_t %>%
                 transform(name = factor(name, levels = dd$labels$label)),
               aes(x = Coord, y = name, col = AF)) +
  geom_tile(linewidth = 1) +
  # scale_color_viridis_c(option = "E") +
  coord_cartesian(xlim = c(23000000, 26000000)) +
  scale_x_continuous(labels = scales::label_comma(scale = 1/1000000)) + # 単位M(1000000)をつけるため値自体は1/1000000する
  scale_color_gradientn(colours = pals::ocean.deep(100), name = "Allele \nfrequency") +
  # xlab("Linkage group 17") +
  xlab("Linkage group 3 (Mb)") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    # axis.ticks = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(),
    panel.border = element_blank(),
    plot.margin = unit(c(0,0,0,0), "null"),
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  )
# g5c1

g5c1_pbs2 <- g5c1_pbs / g5c1 + plot_layout(heights = c(1, 2))
ggsave("../figures/tmp/Figure5c_pbs.pdf", g5c1_pbs2, w = 6, h = 3)

 

library(cowplot)
g5c1_0 <- plot_spacer() / plot_spacer() / g5c1_dend + plot_layout(heights = c(1, 6, 6))
g5c1_1 <- g5c1_gene / g5c1_pbs / g5c1 + plot_layout(heights = c(1, 6, 6))
# g5c <- g5c_0 + g5c_1 + plot_layout(widths = c(1, 5))
g5c <- cowplot::plot_grid(g5c1_0, g5c1_1, ncol = 2, rel_widths = c(1, 6), align = "h")

g5c
ggsave("../figures/tmp/Figure5c.pdf", g5c, w=8, h=4)


# BiocManager::install("ggbio")
library(ggbio)
# BiocManager::install("GenomicRanges")
library(GenomicRanges)
# BiocManager::install("rtracklayer")
library(rtracklayer) # imports gff/bed/wig
gff2 <- import("../../data/analyzed_data/Ref/GCF_020424385.1_ASM2042438v2_genomic.gtf") # contains gene/mRNA/exon/CDS
region <- c(23000000, 26000000)
gr <- GRanges(seqnames = "LG03_CM035475.1", IRanges(region[1], region[2]), strand = "*")
# wh <- genesymbol[c("LOC123758689", "LOC123761000", "LOC123758696")]
# wh <- range(wh, ignore.strand = TRUE)
p <- autoplot(gff2, which=gr) +
  theme_bw()
# ggsave("../../ザリガニ論文_Sato_Makino/figures/tmp/Figure5c_gene.pdf", p, w = 4, h = 2)
pdf("../figures/tmp/Figure5c_gene.pdf", w = 4, h = 2)
print(p)
dev.off()

BiocManager::install("Sushi")
library(Sushi)
bdata <- df_gtf %>%
  filter(X1 == "NC_059573.1",
         X4 > region[1],
         X5 < region[2],
         X3 == "exon") %>%
  dplyr::rename(chrom = X1,
                start = X4,
                stop = X5, 
                gene = geneID,
                score = X8,
                strand = X7,
                type = X3) %>%
  mutate(strand = if_else(strand == "+", "1", "-1")) %>%
  dplyr::select(chrom, start, stop, gene,
                score, strand, type) %>%
  ungroup()

p <- Sushi::plotGenes(bdata, 
                      "NC_059573.1",
                      23000000,
                      26000000,
                      types = bdata$type,
          maxrows = 50,
          bheight = 0.08,
          plotgenetype = "box",
          bentline = FALSE,
          col = "brown",
          labeloffset = .2,
          fontsize = 0.9,
          arrowlength = 0.025,
          labeltext = TRUE)

df_gtf2 <- df_gtf %>%
  

df_gtf_exon <- df_gtf %>%
  filter(X4 > region[1], X5 < region[2], X3 == "exon") %>%
  filter(X1 == "NC_059573.1")
xmin <- df_gtf_exon$X4
xmax <- df_gtf_exon$X5
df_gtf_gene <- df_gtf %>%
  filter(X4 > region[1], X5 < region[2], X3 == "gene") %>%
  filter(X1 == "NC_059573.1")

#### Figure 5d two genes ####
df_tpm <- read.table("../../data/analyzed_data/RNAseq/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104.tsv", header = TRUE) %>%
  dplyr::rename(GeneID = ID) %>%
  pivot_longer(cols = contains("_"), values_to = "TPM") %>%
  separate(col = name, into = c("group", "no"), sep = "_") %>%
  mutate(pop = str_sub(group, 1, 3),
         group = case_when(group == "SEN0" ~ "Sendai Day 0",
                           group == "SEN7" ~ "Sendai Day 7",
                           group == "TON0" ~ "Sapporo Day 0",
                           group == "TON7" ~ "Sapporo Day 7",
                           group == "TON30" ~ "Sapporo Day 31",))

g5d_2genes <- ggplot(df_tpm %>%
                    filter(GeneID %in% c("LOC123756530", "LOC123758570")) %>%
                    transform(group = factor(group, levels = c("Sendai Day 0", "Sendai Day 7", "Sapporo Day 0", "Sapporo Day 7", "Sapporo Day 31"))),
                  aes(x = group, y = log2(TPM+1), col = pop, fill = pop)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               aes(col = pop), width = 0, linewidth = 1.5) +
  scale_x_discrete(label = c("Day 0", "Day 7", "Day 0", "Day 7", "Day 31")) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim=c(0,NA)) +
  scale_color_manual(values = c("#a99e93", "#ebe1a9")) +
  scale_fill_manual(values = c("#a99e93", "#ebe1a9")) +
  ylab("Log2(TPM+1)") +
  facet_wrap( ~ GeneID, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        legend.position = "none",
        strip.background = element_blank())

g5d_2genes
ggsave("../figures/tmp/Figure5d_2genes.pdf", g5d_2genes, width = 2.6, height = 3.4)

