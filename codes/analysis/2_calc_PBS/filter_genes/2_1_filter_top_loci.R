library(tidyverse)
library(ggpubr)
library(iterators)

setwd("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/")

#### 4 < cov < 16 ####
##### 1kbp #####
##Aomori kamakura2 tonden2
df_PBS <- read.table("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/manuscript_tmp/data/PBS/rawdata/pooled_7pops_indelfiltered_mincov4_maxcov16_1kbp_annotatedCHR_Aomori_kamakura2_tonden2_PBS.tsv", h = T, stringsAsFactors = FALSE)

don_PBS <- df_PBS %>%
  
  # filter(df_PBS != "na") %>%
  
  # Compute CHR_tmpomosome size
  group_by(CHR_tmp) %>% 
  summarise(CHR_tmp_len=max(BP)) %>% 
  
  # Calculate cumulative position of each CHR_tmpomosome
  mutate(tot=cumsum(as.numeric(CHR_tmp_len))-CHR_tmp_len) %>%
  select(-CHR_tmp_len) %>%
  
  # Add this info to the initial dataset
  left_join(df_PBS, ., by=c("CHR_tmp"="CHR_tmp")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR_tmp, BP) %>%
  mutate(BPcum=BP+tot)

axisdf <- don_PBS %>% group_by(CHR_tmp) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g21 <- ggplot(don_PBS, aes(x=BPcum, y=as.numeric(PBS_tonden2), color=as.factor(CHR_tmp))) +
  
  # Show all points
  geom_point(size=1.3,alpha=0.8) +
  scale_color_manual(values = rep(c("grey", "black"), 7536 )) +
  #  scale_size_manual(values = c(1.3, 2.4)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  #  scale_fill_viridis_d() +
  #  scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  geom_hline(yintercept = quantile(na.omit(don_PBS$PBS_tonden2),0.99), linewidth = 1, col = "black", linetype = "dotted") +
  geom_hline(yintercept = quantile(na.omit(don_PBS$PBS_tonden2),0.999), linewidth = 1, col = "darkred", linetype = "dotted") +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR_tmp, breaks = axisdf$center ) +
  #  scale_y_continuous(breaks=seq(0,1,by=0.2)) +     # remove space between plot area and x axis
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Linkage Group") +
  ylab("PBS (Tonden)") +
  coord_cartesian(ylim = c(0,NA)) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g21
# ggsave("../figures/tmp/PBS/PBS_Aomori_kamakura2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp.png", g21, w=10, h=3)

df_PBS[order(df_PBS$PBS_tonden2, decreasing = T),] %>%
  filter(PBS_tonden2 > quantile(na.omit(df_PBS$PBS_tonden2),0.999)) %>%
  write.table("../data/PBS/1kbp/PBS_Aomori_kamakura2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp-top0.1.tsv",row.names = F, quote = F, sep = "\t")

df_PBS[order(df_PBS$PBS_tonden2, decreasing = T),] %>%
  filter(PBS_tonden2 > quantile(na.omit(df_PBS$PBS_tonden2),0.99)) %>%
  write.table("../data/PBS/1kbp/PBS_Aomori_kamakura2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp-top1.tsv",row.names = F, quote = F, sep = "\t")


##Aomori okinawa2 tonden2
df_PBS <- read.table("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/manuscript_tmp/data/PBS/rawdata/pooled_7pops_indelfiltered_mincov4_maxcov16_1kbp_annotatedCHR_Aomori_okinawa2_tonden2_PBS.tsv", h = T, stringsAsFactors = FALSE)

don_PBS <- df_PBS %>%
  
  # filter(df_PBS != "na") %>%
  
  # Compute CHR_tmpomosome size
  group_by(CHR_tmp) %>% 
  summarise(CHR_tmp_len=max(BP)) %>% 
  
  # Calculate cumulative position of each CHR_tmpomosome
  mutate(tot=cumsum(as.numeric(CHR_tmp_len))-CHR_tmp_len) %>%
  select(-CHR_tmp_len) %>%
  
  # Add this info to the initial dataset
  left_join(df_PBS, ., by=c("CHR_tmp"="CHR_tmp")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR_tmp, BP) %>%
  mutate(BPcum=BP+tot)

axisdf <- don_PBS %>% group_by(CHR_tmp) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g22 <- ggplot(don_PBS, aes(x=BPcum, y=as.numeric(PBS_tonden2), color=as.factor(CHR_tmp))) +
  
  # Show all points
  geom_point(size=1.3,alpha=0.8) +
  scale_color_manual(values = rep(c("grey", "black"), 7536 )) +
  #  scale_size_manual(values = c(1.3, 2.4)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  #  scale_fill_viridis_d() +
  #  scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  geom_hline(yintercept = quantile(na.omit(don_PBS$PBS_tonden2),0.99), linewidth = 1, col = "black", linetype = "dotted") +
  geom_hline(yintercept = quantile(na.omit(don_PBS$PBS_tonden2),0.999), linewidth = 1, col = "darkred", linetype = "dotted") +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR_tmp, breaks = axisdf$center ) +
  #  scale_y_continuous(breaks=seq(0,1,by=0.2)) +     # remove space between plot area and x axis
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Linkage Group") +
  ylab("PBS (Tonden)") +
  coord_cartesian(ylim = c(0,NA)) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g22
# ggsave("../figures/tmp/PBS/PBS_Aomori_okinawa2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp.png", g22, w=10, h=3)

df_PBS[order(df_PBS$PBS_tonden2, decreasing = T),] %>%
  filter(PBS_tonden2 > quantile(na.omit(df_PBS$PBS_tonden2),0.999)) %>%
  write.table("../data/PBS/1kbp/PBS_Aomori_okinawa2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp-top0.1.tsv",row.names = F, quote = F, sep = "\t")

df_PBS[order(df_PBS$PBS_tonden2, decreasing = T),] %>%
  filter(PBS_tonden2 > quantile(na.omit(df_PBS$PBS_tonden2),0.99)) %>%
  write.table("../data/PBS/1kbp/PBS_Aomori_okinawa2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp-top1.tsv",row.names = F, quote = F, sep = "\t")


##kamakura2 okinawa2 tonden2
df_PBS <- read.table("/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/manuscript_tmp/data/PBS/rawdata/pooled_7pops_indelfiltered_mincov4_maxcov16_1kbp_annotatedCHR_kamakura2_okinawa2_tonden2_PBS.tsv", h = T, stringsAsFactors = FALSE)

don_PBS <- df_PBS %>%
  
  # filter(df_PBS != "na") %>%
  
  # Compute CHR_tmpomosome size
  group_by(CHR_tmp) %>% 
  summarise(CHR_tmp_len=max(BP)) %>% 
  
  # Calculate cumulative position of each CHR_tmpomosome
  mutate(tot=cumsum(as.numeric(CHR_tmp_len))-CHR_tmp_len) %>%
  select(-CHR_tmp_len) %>%
  
  # Add this info to the initial dataset
  left_join(df_PBS, ., by=c("CHR_tmp"="CHR_tmp")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR_tmp, BP) %>%
  mutate(BPcum=BP+tot)

axisdf <- don_PBS %>% group_by(CHR_tmp) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g23 <- ggplot(don_PBS, aes(x=BPcum, y=as.numeric(PBS_tonden2), color=as.factor(CHR_tmp))) +
  
  # Show all points
  geom_point(size=1.3,alpha=0.8) +
  scale_color_manual(values = rep(c("grey", "black"), 7536 )) +
  #  scale_size_manual(values = c(1.3, 2.4)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  #  scale_fill_viridis_d() +
  #  scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  geom_hline(yintercept = quantile(na.omit(don_PBS$PBS_tonden2),0.99), linewidth = 1, col = "black", linetype = "dotted") +
  geom_hline(yintercept = quantile(na.omit(don_PBS$PBS_tonden2),0.999), linewidth = 1, col = "darkred", linetype = "dotted") +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR_tmp, breaks = axisdf$center ) +
  #  scale_y_continuous(breaks=seq(0,1,by=0.2)) +     # remove space between plot area and x axis
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Linkage Group") +
  ylab("PBS (Tonden)") +
  coord_cartesian(ylim = c(0,NA)) +
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g23
# ggsave("../figures/tmp/PBS/PBS_kamakura2_okinawa2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp.png", g23, w=10, h=3)

df_PBS[order(df_PBS$PBS_tonden2, decreasing = T),] %>%
  filter(PBS_tonden2 > quantile(na.omit(df_PBS$PBS_tonden2),0.999)) %>%
  write.table("../data/PBS/1kbp/PBS_kamakura2_okinawa2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp-top0.1.tsv",row.names = F, quote = F, sep = "\t")

df_PBS[order(df_PBS$PBS_tonden2, decreasing = T),] %>%
  filter(PBS_tonden2 > quantile(na.omit(df_PBS$PBS_tonden2),0.99)) %>%
  write.table("../data/PBS/1kbp/PBS_kamakura2_okinawa2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp-top1.tsv",row.names = F, quote = F, sep = "\t")


