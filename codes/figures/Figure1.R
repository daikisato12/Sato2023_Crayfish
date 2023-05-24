#### load libraries ####
library(tidyverse)
library(data.table) #fread
library(rlist)
library(pals)
library(ggpubr)
library(maps)
library(scales)
library(igraph)
library(ggnetwork)
# install.packages("ggmap")
# library(ggmap)
# install.packages("countrycode")
library(countrycode)
# devtools::install_github("uribo/jpndistrict")
library(jpndistrict)
# library(iterators)
# library(gtools) #mixedsort

#### Figure 1a Map ####
world <- ggplot2::map_data("world")
g1_US <- world %>%
  filter(region == "USA") %>% 
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightgray", colour = "black", linewidth = 0.1) +
  coord_cartesian(xlim = c(-130, -60), ylim = c(20, 50)) +
  xlab("Longitude (°)") +
  ylab("Latitude (°)") +
  theme_minimal()
# g1_US
ggsave("../figures/tmp/Figure1a_US.pdf", g1_US, w = 5, h = 3)

g1_JP <- world %>%
  filter(region == "Japan") %>% 
  ggplot(aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightgray", colour = "black", linewidth = 0.1) +
  # coord_cartesian(xlim = c(-130, -60), ylim = c(20, 50)) +
  xlab("Longitude (°)") +
  ylab("Latitude (°)") +
  theme_minimal()
# g1_JP
ggsave("../figures/tmp/Figure1a_JP.pdf", g1_JP, w = 3, h = 4)

g1a <- ggarrange(g1_JP, g1_US, ncol = 2, widths = c(1, 2.8), align = "v")
g1a
ggsave("../figures/tmp/Figure1a_map.pdf", g1a, w = 10, h = 3)


#### Figure 1b PSMC ####
"+" = function(e1, e2){ #make function just to make it easy to connect strings
  if(is.character(c(e1, e2))){
    paste(e1, e2, sep="")
  }else{
    base::"+"(e1, e2)
  }
}

#files <- list.files("./PSMC/") %>% str_subset("txt")
files <- list.files("../../data/analyzed_data/PSMC/Ref_Pcla_dep8sites/bs/", recursive=T) %>% 
  str_subset(".results.bs.0.txt")
name_list <- c("atchafalaya" = "atchafalaya",
               "sendai" = "sendai",
               "SRR14457223" = "china1",
               "SRR14457234" = "china2",
               "SRR14457235" = "china3",
               "SRR5115141" = "virginalis",
               "SRR5115151" = "fallax",
               "SRR5115153" = "alleni",
               "zonangulus" = "zonangulus")
df_psmc <- data.frame()
for (file in files) {
  pop <- str_split(file,"\\.") %>% map_chr(1)
  pop <- name_list[pop]
  
  testdf <- read.delim("../../data/analyzed_data/PSMC/Ref_Pcla_dep8sites/bs/" + file, h=F)
  YA <- testdf$V1
  Ne <- testdf$V2
  n.points <- nrow(testdf)
  
  YearsAgo<-c(as.numeric(rbind(YA[-n.points],YA[-1])), YA[n.points])
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])), Ne[n.points])
  
  df_t <- bind_cols(YearsAgo, Ne) %>%
    dplyr::rename("YearsAgo" = ...1, "Ne" = ...2) %>%
    mutate(Population = pop)
  
  df_psmc <- bind_rows(df_psmc, df_t)
  
}

spe_label <- c(expression(italic("P. alleni")), 
               expression(paste(italic("P. clarkii"), " (Atchafalaya)")),
               expression(paste(italic("P. clarkii"), " (China 1)")),
               expression(paste(italic("P. clarkii"), " (China 2)")),
               # expression(paste(italic("P. clarkii"), " (Sendai)")),
               expression(italic("P. fallax")),
               expression(italic("P. virginalis")),
               expression(italic("P. zonangulus"))) 

g1b <- ggplot(df_psmc %>%
                mutate(Population = factor(Population,
                                           levels = c("alleni", "atchafalaya", "china1", "china2", "china3", "sendai", "fallax", "virginalis", "zonangulus"))) %>%
                filter(!Population %in% c("sendai", "china1")), 
              aes(x = YearsAgo, y = Ne)) +
  annotate("rect", xmin = 1.17e4, xmax = 1.15e5, ymin = -Inf, ymax = Inf,  fill = "lightblue", alpha=.6) +
  geom_line(aes(x = YearsAgo, y = Ne, color = Population)) +
  coord_cartesian(xlim = c(5e2, 2e6), ylim=c(0, 40)) + # 40
  annotation_logticks(sides = "b", long = unit(2,"mm"), mid = unit(1,"mm"), size = .2) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), breaks = 10^(3:7)) +
  scale_color_manual(values = c("#b28c6e", "#d9333f", "#a25768", "#887f7a", "#674598", "#0094c8", "#00a381", "#b8d200"), labels = spe_label) + #tol.rainbow(8)
  # labs(color = "Species") +
  xlab(expression(paste("Years ago (g = 1, µ = 4.59 × ", 10^{-9}, ")"))) +
  ylab(expression(paste("Effective population size (×", 10^{4}, ")"))) +
  theme_bw() +
  theme(legend.text.align = 0,
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.54, 0.8)) +
  guides(color = guide_legend(ncol = 2))
g1b
ggsave("../figures/tmp/Figure1b_PSMC_dep8.pdf", g1b, width=4, height=3)


#### Figure 1c Pi 5kbp ####
list_files <- list.files(paste0("../../data/analyzed_data/Pi/"), pattern = "*.subsampled.pi.5kb", full.names = TRUE, recursive=T)
df_pi <- data.frame()
for (i in 1:length(list_files)){
  # i <- 1
  sample_name <- list_files[i] %>%
    basename() %>%
    str_split("\\.") %>%
    map_chr(1)
  df_pi <- bind_rows(df_pi, 
                     data.frame(Sample = sample_name,
                                read.table(list_files[i], header = F))
  )
}
colnames(df_pi) <- c("Sample", "Chr", "Coord", "N_SNPs", "Cover_fraction", "pi")

df_pi_2 <- df_pi %>%
  filter(pi != "na") %>%
  mutate(pi = as.numeric(pi)) %>%#,
  mutate(Sample = case_when(Sample == "NewOrleans" ~ "New Orleans",
                            Sample == "triunfo2" ~ "Triunfo",
                            Sample == "atchafalaya2" ~ "Atchafalaya",
                            Sample == "kamakura2" ~ "Kamakura",
                            Sample == "tonden2" ~ "Sapporo",
                            Sample == "Aomori" ~ "Aomori",
                            Sample == "okinawa2" ~ "Okinawa"),
         Country = case_when(Sample == "Kamakura" ~ "Japan",
                             Sample == "Sapporo" ~ "Japan",
                             Sample == "Aomori" ~ "Japan",
                             Sample == "Okinawa" ~ "Japan",
                             TRUE ~ "US")) %>%
  transform(Sample = factor(Sample, levels = c("New Orleans", "Triunfo", "Atchafalaya", "Kamakura", "Sapporo", "Aomori", "Okinawa")))

g1c <- ggpubr::ggerrorplot(df_pi_2,
              x = "Sample", y = "pi", color = "Country", 
              desc_stat = "mean_ci", size = 1, alpha = 0.5) +
  scale_y_continuous(labels = scales::label_number(scale = 1000),
  breaks = seq(0,0.01,0.0002)) +
  scale_color_manual(values = c("#e9546b", "#84a2d4")) +
  ylab(expression(italic(pi) (10^-3))) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.title.x = element_blank(),
        legend.position = "none")
g1c
ggsave("../figures/tmp/Figure1c_Pi_5kbp.pdf", g1c, w = 3, h = 3)


##### analysis #####
pairwise.wilcox.test(g = df_pi_2$Sample, x = df_pi_2$pi)

#### Figure 1d Fst ####
combinations_list <- c("Aomori-Atchafalaya", "Aomori-Kamakura", "Aomori-New Orleans", "Aomori-Okinawa", "Aomori-Tonden", "Aomori-Triunfo",
                       "Atchafalaya-Kamakura", "Atchafalaya-New Orleans", "Atchafalaya-Okinawa", "Atchafalaya-Tonden", "Atchafalaya-Triunfo",
                       "Kamakura-New Orleans", "Kamakura-Okinawa", "Kamakura-Tonden", "Kamakura-Triunfo", 
                       "New Orleans-Okinawa", "New Orleans-Tonden", "New Orleans-Triunfo", 
                       "Okinawa-Tonden", "Okinawa-Triunfo",
                       "Tonden-Triunfo")


df_fst <- read.table("../data/Fst/pooled_7pops_indelfiltered_mincov4_maxcov16_5kbp.fst", header = F) %>%
  `colnames<-`(c("Chr", "Coord", "N_SNPs", "Fraction", "Coverage", combinations_list)) %>%
  pivot_longer(cols = contains("-"), names_to = "Pop1-Pop2", values_to = "Fst") %>% 
  mutate(Fst = str_sub(Fst, 5) %>% as.numeric()) %>%
  group_by(`Pop1-Pop2`) %>%
  summarize(Fst = mean(Fst, na.rm = T)) %>%
  separate(col = `Pop1-Pop2`, into = c("Pop1", "Pop2"), sep = "-") #%>%

##### Heatmap visualization #####
g1d_heat <- ggplot(df_fst, aes(x = Pop1, y = Pop2)) +
  geom_tile(aes(fill = Fst)) +
  scale_fill_gradientn(colors = rev(pals::brewer.brbg(100)), 
                        name = expression(italic(F)[ST])) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        axis.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.position = c(1, 0), 
        legend.justification = c(1, -0.1))

ggsave("../figures/tmp/Figure1d_Fst_heatmap.pdf", g1d_heat, w = 3, h = 3)


##### Network visualization #####
# install.packages("ggnewscale")
library(ggnewscale)

x <- list("Atchafalaya" = 0, "New Orleans" = -0.6, "Triunfo" = 0, "Kamakura" = 1.2, "Aomori" = 1.8, "Tonden" = 1.8, "Okinawa" = 1.2)
y <- list("Atchafalaya" = 0, "New Orleans" = 0.6, "Triunfo" = 1.2, "Kamakura" = 0, "Aomori" = 0.4, "Tonden" = 0.8, "Okinawa" = 1.2)
df_fst_network <- df_fst %>%
  bind_cols(data.frame(x = unlist(x[df_fst$Pop1]),
                       xend = unlist(x[df_fst$Pop2]),
                       y = unlist(y[df_fst$Pop1]),
                       yend = unlist(y[df_fst$Pop2]))) %>%
  mutate(Country = if_else(Pop1 %in% c("Atchafalaya", "New Orleans", "Triunfo"), "U.S.", "Japan"))


g1d_net <- ggplot(df_fst_network, aes(x = x, y = y, xend = xend, yend = yend)) +
  ggnetwork::geom_edges(aes(linewidth = Fst, color = Fst)) +
  scale_color_gradientn(colors = rev(pals::brewer.brbg(100)), 
                        name = expression(italic(F)[ST])) +
  new_scale_colour() +
  ggnetwork::geom_nodes(aes(color = Country), size = 10) +
  scale_color_manual(values = c("#e9546b", "#84a2d4")) +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "top")

ggsave("../figures/tmp/Figure1d_Fst_network.pdf", g1d_net, w = 5.6, h = 6)



#### Figure 1e ####
# Created by iTOL (https://itol.embl.de)
