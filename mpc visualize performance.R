library(ggplot2)
library(ggpubr)
library(dplyr)
library(viridis)

data_ch = read.csv('competitive_hierarchy_mpc_stats_7n.csv') %>% mutate(Hypothesis='Competitive hierarchy')
data_ls = read.csv('limiting_similarity_mpc_stats_7n.csv') %>% mutate(Hypothesis='Limiting similarity')

data_all = rbind(data_ch, data_ls) %>% select(-X)

g_perf = ggplot(data_all, aes(y=Controller_Success,
                     color=Target_Trait_Composition, 
                     x=Initial_Richness,
                     group=Target_Trait_Composition,
                     ymin=Controller_Success-Confidence_Interval_95/2,
                     ymax=Controller_Success+Confidence_Interval_95/2)) +
  geom_point() +
  geom_line() + 
  geom_errorbar(width=0.2) +
  facet_wrap(~Hypothesis) +
  ylim(0,1) +
  theme_bw() +
  scale_color_viridis(name='Community-weighted\nmean trait value',lim=c(0,1)) +
  ylab("Controller success rate") +
  xlab("Initial species richness")
  


ggsave(g_perf,file='g_mpc_success_rate.png',width=7,height=4)
