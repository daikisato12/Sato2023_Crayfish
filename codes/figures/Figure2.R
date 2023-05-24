#### load libraries ####
library(tidyverse)
library(ggpubr)
# install.packages("survival")
library(survival)
# install.packages("rstatix")
# library(rstatix)
# install.packages("survminer")
library(survminer)
# install.packages("tsibble")
library(tsibble)

#### Figure 2a Temperature ####
df_temp <- read.csv("../../data/analyzed_data/Temperature/data_rev.csv") %>%
  separate(Month, sep = "-", into = c("month", "year")) %>%
  mutate(year = if_else(parse_number(year) < 21, paste0("20", year), paste0("19", year)),
         yearmonth = tsibble::yearmonth(paste0(year, " ", month))) %>%
  pivot_longer(cols = contains("_"), names_sep = "_", names_to = c("pop", "var", "unit")) %>%
  arrange(yearmonth, pop) %>%
  transform(month = factor(month, levels = c("Jan", "Feb", "Mar", "Apr",
                                             "May", "Jun", "Jul", "Aug",
                                             "Sep", "Oct", "Nov", "Dec")),
            pop = factor(pop, levels = c("Sendai", "Sapporo")),
            var = factor(var, levels = c("MaximumTemperature", "MinimumTemperature",
                                         "AmountSnowfall", "DaysSnowfall")))

g2a <- ggplot(df_temp %>%
                filter(var %in% c("AmountSnowfall", "MinimumTemperature")), 
              aes(x = month, y = value, col = pop, group = pop)) +
  stat_summary(fun.data = mean_sd, na.rm =TRUE, 
               geom = "errorbar", width = 0, linewidth = .4) +
  stat_summary(fun = mean, na.rm = TRUE,
               geom = "line", 
               size = .8) +
  ylab("") +
  scale_color_manual(values = c("#a99e93", "#ebe1a9"), labels = c("Sendai", "Sapporo")) +
  facet_wrap(~ var, nrow = 4, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.title = element_blank())
g2a
ggsave("../figures/tmp/Figure2a_temperature.pdf", g2a, w = 5, h = 3)


#### Figure 2b survival analysis ####

df_surv <- read.delim("../../data/analyzed_data/SurvivalExperiment/survival.txt")
fit <- survfit(Surv(survival, censoring) ~ pop, data = df_surv, type = "kaplan-meier") #Kaplan-Meier curve

g2b <- ggsurvplot(
  fit,                     
  data = df_surv,             
  # risk.table = TRUE,       
  pval = TRUE, 
  pval.method = TRUE,
  # conf.int = TRUE,
  xlab = "Days", 
  xlim = c(0, 35),
  ggtheme = theme_bw(), 
  # risk.table.y.text.col = T, 
  # risk.table.y.text = FALSE, 
  legend.labs = c("Sendai (n = 74)", "Tonden (n = 29)"),
  legend.title = element_blank(),
  palette = c("#a99e93", "#ebe1a9")

  # censor.shape = 16,
)

g2b
ggsave("../figures/tmp/Figure2b_survival.pdf", width=3.5, height=4)
