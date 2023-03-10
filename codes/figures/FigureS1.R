#### load libraries ####
library(tidyverse)
setwd("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/codes/")

#### Figure S1 TPM ####

df_table1 <- read.table("../figures/tmp/Table1.tsv", header = TRUE, sep = "\t") %>%
  # tibble::rownames_to_column() %>%
  separate(col = Crayfish.gene.ID, into = paste0("LOC", seq(1,16))) %>%
  pivot_longer(cols = contains("LOC"), values_to = "GeneID") %>%
  na.omit() %>%
  dplyr::select(!name) %>%
  filter(Fly.ortholog %in% c("Ank2", "CG4115", "CG6055", "Col4a1", "Cp1",
                             "dpy", "Meics", "S", "Spn38F", "Spn88Ea", "Tep2"))
# dplyr::rename(Genename = FlyID,
# FlyID = rowname)

df_tpm <- read.table("../data/RNA-seq/unique_mapping/rawdata/crayfish_hisat2_stringtie_uniq_tpm_genes_2n1tpm20230104.tsv", header = TRUE) %>%
  dplyr::rename(GeneID = ID) %>%
  left_join(df_table1) %>%
  pivot_longer(cols = contains("_"), values_to = "TPM") %>%
  separate(col = name, into = c("group", "no"), sep = "_") %>%
  mutate(pop = str_sub(group, 1, 3),
         group = case_when(group == "SEN0" ~ "Sendai Day 0",
                           group == "SEN7" ~ "Sendai Day 7",
                           group == "TON0" ~ "Sapporo Day 0",
                           group == "TON7" ~ "Sapporo Day 7",
                           group == "TON30" ~ "Sapporo Day 31",))

# g3d_Spn <- ggplot(df_tpm %>%
#                     filter(str_detect(Genename, "Spn")) %>%
#                     transform(group = factor(group, levels = c("Sendai Day 0", "Sendai Day 7", "Sapporo Day 0", "Sapporo Day 7", "Sapporo Day 31"))), 
#                   aes(x = group, y = log2(TPM+1), col = pop, fill = pop)) +
#   geom_bar(stat = "summary", fun = "mean", width = 0.6) +
#   stat_summary(fun.data = "mean_se", geom = "errorbar", 
#                aes(col = pop), width = 0, linewidth = 1.5) +
#   scale_x_discrete(label = c("Day 0", "Day 7", "Day 0", "Day 7", "Day 31")) +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_cartesian(ylim=c(0,NA)) +
#   scale_color_manual(values = c("#a99e93", "#ebe1a9")) +
#   scale_fill_manual(values = c("#a99e93", "#ebe1a9")) +
#   ylab("Log2(TPM+1)") +
#   facet_wrap(Genename ~ GeneID, scales = "free_y") +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
#         legend.position = "none",
#         strip.background = element_blank())
# 
# g3d_Spn
# ggsave("../figures/tmp/Figure3d_Spn.pdf", g3d_Spn, width = 6, height = 6)


df_PBS <- read.table("../data/PBS/HAF_PBS_5kbp_step1kbp-top0.1_2combsharedloci_genelist.txt", header = FALSE) %>%
  rename(GeneID = V1) %>%
  mutate(selection = "selected")

df_S1 <- df_tpm %>%
  filter(!is.na(Fly.ortholog)) %>%
  group_by(Fly.ortholog) %>%
  mutate(gene_no = as.integer((row_number()-1)/15)) %>%
  ungroup() %>%
  left_join(df_PBS) %>%
  mutate(selection = if_else(is.na(selection), "not selected", selection)) %>%
  transform(group = factor(group, levels = c("Sendai Day 0", "Sendai Day 7", "Sapporo Day 0", "Sapporo Day 7", "Sapporo Day 31")))

gS1 <- ggplot(df_S1, aes(x = group, y = log2(TPM+1), col = pop, fill = pop)) +
  geom_rect(aes(fill = selection, alpha = selection), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", 
               aes(col = pop), width = 0, linewidth = 1) +
  # scale_x_discrete(label = c("Day 0", "Day 7", "Day 0", "Day 7", "Day 31")) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim=c(0,NA)) +
  scale_color_manual(values = c("#a99e93", "#ebe1a9")) +
  scale_fill_manual(values = c("selected" = "#eebbcb", "not selected" = NA, "SEN" = "#a99e93", "TON" = "#ebe1a9")) +
  scale_alpha_manual(values = c("selected" = 0.2, "not selected" = 0)) +
  ylab(expression(Log[2](TPM+1))) +
  facet_grid(Fly.ortholog ~ gene_no, scales = "free_y") +
  geom_text(data = df_S1 %>%
              group_by(Fly.ortholog, gene_no, GeneID) %>%
              slice(1),
            mapping = aes(x = -Inf, y = Inf, label = GeneID), hjust = -0.1, vjust = 2, color = "black", size = 2) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

gS1
ggsave("../figures/tmp/FigureS1.pdf", gS1, width = 14, height = 8)

# col_order <- c(paste0("Sendai Day 0 rep", seq(1,3)),
#                paste0("Sendai Day 7 rep", seq(1,3)),
#                paste0("Sapporo Day 0 rep", seq(1,3)),
#                paste0("Sapporo Day 7 rep", seq(1,3)),
#                paste0("Sapporo Day 31 rep", seq(1,3)))

# modify_gt3 <- function(x){
#   x <- ifelse(x > 3, 3, x)
#   return(x)
# }
# 
# df_tpm_dup <- df_tpm %>%
#   filter(!is.na(Genename)) %>%
#   mutate(sample = paste0(group, " rep", no),
#          rowname = paste0(Genename, "_", GeneID)) %>%
#   group_by(GeneID, Genename, rowname, sample) %>%
#   summarize(TPM = mean(TPM, na.rm = T)) %>%
#   ungroup() %>%
#   dplyr::select(rowname, sample, TPM) %>%
#   pivot_wider(names_from = sample, values_from = TPM) %>%
#   arrange(rowname) %>%
#   tibble::column_to_rownames()
# 
# df_tpm_dup <- scale(log2(df_tpm_dup+1))
# df_tpm_dup <- df_tpm_dup[,col_order] %>%
#   modify_gt3()
# 
# # df_tpm_dup <- df_tpm_dup[c("GeneID", "Genename", col_order)]
# 
# # library(heatmap3)
# colors <- colByValue(df_tpm_dup, 
#                      col = colorRampPalette(c("#a2d7dd", "#d6c6af", "#c97586"))(51),
#                      breaks = c(seq(-3,3,length=52)))
# dynamicColors <- as.data.frame(colors)
# rownames(dynamicColors) <- rownames(df_tpm_dup)
# colnames(dynamicColors) <- colnames(df_tpm_dup)
# dynamicColors[is.na(dynamicColors)] <- "#ffffff"
# 
# df_fig3d <- dynamicColors %>%
#   tibble::rownames_to_column() %>%
#   pivot_longer(cols = !rowname, names_to = "var", values_to = "color") %>%
#   transform(rowname = factor(rowname, levels = rev(rownames(df_tpm_dup))),
#             var = factor(var, levels = colnames(df_tpm_dup)))
# 
# g_fig3d <- ggplot(df_fig3d, aes(x = var, y = rowname, fill = color)) +
#   geom_tile() +
#   scale_fill_identity() +
#   scale_x_discrete(position = "top") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1, size = 6),
#         axis.text.y = element_text(size = 6),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(), 
#         panel.grid = element_blank())
# g_fig3d

