library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(viridis)
library(ggrepel)

# add in time series
data_statevars = read.csv("preds_statevars.csv",header=FALSE)
names(data_statevars) = LETTERS[1:ncol(data_statevars)]
data_perturbations = read.csv("preds_perturbations.csv",header=FALSE)
names(data_perturbations) = paste("u",LETTERS[1:ncol(data_statevars)],sep="")

data_time = read.csv("preds_time.csv",header=FALSE)
names(data_time) = "t"
data_aux = read.csv("preds_aux.csv",header=FALSE)
names(data_aux) = c("Blank","T_cwm","N_tot","richness")

data_traits = read.csv("preds_traits.csv",header=F) %>% rename(trait.value=V1) %>% mutate(species=LETTERS[1:ncol(data_statevars)])
data_traits = rbind(data_traits, data_traits %>% mutate(species=paste("u",species,sep="")))

data_ts = cbind(data_time, data_statevars, data_perturbations, data_aux)

data_ts_melted = melt(data_ts,id.vars=c("t")) %>%
  left_join(data_traits %>% rename(variable=species))

vars_perturbed = data_ts_melted %>% group_by(variable) %>% summarize(mean=mean(value)) %>% filter(variable %in% paste("u",LETTERS[1:ncol(data_statevars)],sep="")) %>% filter(mean!=0)



g_perturbation = ggplot(data=data_ts_melted %>% filter(variable %in% paste("u",LETTERS[1:ncol(data_statevars)],sep="") & variable %in% vars_perturbed$variable)) +
  geom_line(aes(x=t,y=value,col=trait.value,group=variable),size=0.5) + 
  theme_bw() +
  xlab("Time (t)") + ylab("Perturbation (u)") +
  geom_label_repel(data=data_ts_melted %>% filter(t==max(t) & variable %in% paste("u",LETTERS[1:ncol(data_statevars)],sep="") & variable %in% vars_perturbed$variable),aes(x=t,y=value,label=variable,color=trait.value),max.overlaps=Inf) +
  ggtitle('Controlled variable') +
  scale_color_viridis(name='Trait value',lim=c(0,1)) +
  scale_x_sqrt()

g_controlled = ggplot(data=data_ts_melted %>% filter(variable %in% LETTERS[1:ncol(data_statevars)])) +
  geom_line(aes(x=t,y=value,col=trait.value,group=variable),size=1) + 
  theme_bw() +
  xlab("Time (t)") + ylab("Abundance (N)") +
  geom_label_repel(data=data_ts_melted %>% filter(t==max(t) & variable %in% LETTERS[1:ncol(data_statevars)]),aes(x=t,y=value,label=variable,color=trait.value),max.overlaps=Inf) +
  ggtitle('Controlled dynamics') +
  scale_color_viridis(name='Trait value',lim=c(0,1)) +
  scale_x_sqrt()



g_trait = ggplot(data=data_ts_melted %>% filter(variable=="T_cwm")) + 
  geom_line(aes(x=t,y=value),size=1) +
  theme_bw() + 
  xlab("Time (t)") + ylab("Community-weighted mean trait (T)") +
  ggtitle('Trait response') +
  ylim(0,1) +
  geom_hline(yintercept = 0.5,color='red',linetype='dashed') +
  scale_x_sqrt()

g_richness = ggplot(data=data_ts_melted %>% filter(variable=="richness")) + 
  geom_line(aes(x=t,y=value),size=1) +
  theme_bw() + 
  xlab("Time (t)") + ylab("Species richness") +
  ggtitle('Richness response') +
  ylim(0,10) +
  geom_hline(yintercept = 10,color='red',linetype='dashed') +
  scale_x_sqrt() +
  scale_y_continuous(breaks=0:10,limits=c(0,10))

g_all = ggarrange(g_perturbation, g_controlled, g_trait, g_richness, labels='auto', nrow=1,ncol=4,common.legend=TRUE,legend='bottom')

ggsave(g_all,file='g_control_continuous.png',width=9,height=4)
