library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)


data = read.csv("control_solution.csv",header=FALSE)
data_nocontrol = read.csv("no_control_solution.csv",header=FALSE)

names(data) = c("N1","N2","N3","u","t")
names(data_nocontrol) = c("N2","N3","t")

data_melted = melt(data,id.vars=c("t")) %>%
  rename(Species=variable) %>%
  mutate(variable_type=ifelse(Species %in% c("N2","N3"),"Species",ifelse(Species %in% c("N1"),"Species","Control action")))

data_melted_nocontrol = melt(data_nocontrol,id.vars=c("t")) %>%
  rename(Species=variable) %>%
  mutate(variable_type="Species")

g_oc1 = ggplot() +
  geom_line(data=data_melted %>% dplyr::filter(variable_type=="Species"), aes(x=t,y=value,col=Species),size=1) + 
  geom_line(data=data_melted_nocontrol %>% dplyr::filter(variable_type=="Species"), aes(x=t,y=value,col=Species),linetype='dashed') + 
  theme_bw() +
  xlab("Time (t)") +
  ylab("Abundance") +
  ylim(0,40) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73"))

g_oc2 = ggplot(data_melted %>% dplyr::filter(variable_type!="Species"), aes(x=t,y=value,col=Species)) + 
  geom_line(col='black',size=1) +
  theme_bw() +
  xlab("Time (t)") +
  ylab("Perturbation (u)")


wd = getwd()
setwd("../")
source('predict_environmental_change_response.R')
source('simulate_environmental_change.R')
require('deSolve')
require('viridis')
require('ggpubr')

cyclic.change <- simulate_environmental_change(type = "cyclic", cycle_value = 10, cycle_period = 0.08, 
                                               plot.env.change = F, 
                                               n.time.step=200)

cyclic_change_response <- predict_environmental_change_response(nsp = 10,
                                                                env.change = cyclic.change,
                                                                trait.distribution = "uniform",
                                                                mechanism = "niche difference",
                                                                Nmin = 0.01,
                                                                initial.abundance = 0.01,
                                                                growth.rate = "trait",
                                                                extinction = T,
                                                                plot.species.dynamics = T,
                                                                plot.response.diagram = T)

setwd(wd)

g0 = ggplot(cyclic.change, aes(x=time,y=env.dynamic)) +
  geom_line(color='#004000',size=1) + theme_bw() +
  xlab("Time (t)") + ylab("Environment")

g1 = ggplot(cyclic_change_response[[1]], aes(y = n, x = time)) +
  geom_line(aes(col = trait.i, group = i),size=1) +
  scale_color_viridis("Trait value") +
  theme_bw()+
  theme(legend.position = "none") +
  #scale_y_sqrt() +
  xlab("Time (t)") + ylab("Abundance (N)")

g2 = ggplot(cyclic_change_response$detJ, aes(x=time,y=abs(detJ))) +
  geom_line(color='black',size=1) + theme_bw() + ylab("abs(det J)") + xlab("Time (t)")

g_all = ggarrange(
  ggarrange(g0, g1, g2,labels=c("A","B","C"),nrow=1,ncol=3,common.legend=TRUE),
  ggarrange(g_oc1, g_oc2,labels=c("D","E"),nrow=1,ncol=2,common.legend=TRUE)
,nrow=2,ncol=1)

ggsave(g_all,file='fig_8.pdf',width=6,height=5)
ggsave(g_all,file='fig_8.png',width=6,height=5)
