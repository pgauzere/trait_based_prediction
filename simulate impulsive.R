library(ggrepel)

source("predict_demographic_model_parameters.R")
source('simulate_environmental_change.R')
source('predict_environmental_change_response.R')









# do a step change in the environment as the perturbation
tmax = 100
nsp = 10
env.change = data.frame(time=seq(1,tmax,by=0.1),env.dynamic=8)

env.change.perturbed = env.change
env.change.perturbed$env.dynamic[env.change.perturbed$time > 20 & env.change.perturbed$time < 40] = 1

set.seed(1)
trait_values = sort(runif(n=nsp))

set.seed(1)
env_change_response <- predict_environmental_change_response(nsp = nsp,
                                                             env.change = env.change,
                                                             trait.distribution = 'forced',
                                                             trait_values = trait_values,
                                                             mechanism = "niche difference non normalized",
                                                             Nmin = 0.001,
                                                             initial.abundance = 0.01,
                                                             growth.rate = 0.5,
                                                             extinction = T,
                                                             plot.species.dynamics = F,
                                                             plot.response.diagram = F)
set.seed(1)
env_change_response.perturbed <- predict_environmental_change_response(nsp = nsp,
                                                                       env.change = env.change.perturbed,
                                                                       trait.distribution = 'forced',
                                                                       trait_values = trait_values,
                                                                       mechanism = "niche difference non normalized",
                                                                       Nmin = 0.001,
                                                                       initial.abundance = 0.01,
                                                                       growth.rate = 0.5,
                                                                       extinction = T,
                                                                       plot.species.dynamics = F,
                                                                       plot.response.diagram = F)

ecr = env_change_response$comm_dynamic %>% 
  mutate(label=LETTERS[factor(trait.i,levels=sort(unique(trait.i)))])

ecr.perturbed = env_change_response.perturbed$comm_dynamic %>% 
  mutate(label=LETTERS[factor(trait.i,levels=sort(unique(trait.i)))])

g1 = ggplot(ecr, aes(y = n, x = time)) +
  geom_line(aes(col = trait.i, group = i),size=1) +
  scale_color_viridis("Trait value") +
  theme_bw()+
  theme(legend.position = "none") +
  #scale_y_sqrt() +
  xlab("Time (t)") + ylab("Abundance (N)") +
  geom_label_repel(data=ecr %>% filter(time==tmax),aes(x=time,y=n,label=label,color=trait.i),max.overlaps=Inf) +
  ggtitle('Uncontrolled dynamics')


g1.perturbed = ggplot(ecr.perturbed, aes(y = n, x = time)) +
  geom_line(aes(col = trait.i, group = i),size=1) +
  scale_color_viridis("Trait value") +
  theme_bw()+
  theme(legend.position = "none") +
  #scale_y_sqrt() +
  xlab("Time (t)") + ylab("Abundance (N)") +
  geom_label_repel(data=ecr.perturbed %>% filter(time==tmax),aes(x=time,y=n,label=label,color=trait.i),max.overlaps=Inf) +
  ggtitle('Controlled dynamics')

g.env = ggplot(env.change.perturbed,aes(x=time,y=env.dynamic)) +
  geom_line(size=1) +
  geom_line(data=env.change,linetype='dotted',color='gray',size=1) +
  ylim(0,10) +
  theme_bw() +
  xlab("Time (t)") + ylab("Environment (E)") +
  ggtitle('Controlled variable')

g_final = ggarrange(g.env, g1, g1.perturbed,align='hv',common.legend = TRUE,labels='auto',nrow=1,ncol=3,legend = 'bottom')
ggsave(g_final, file='g_control_environment.png',width=9,height=4)














# TRY ANOTHER CASE WHERE WE ADD REMOVE SPECIES
nsp = 10
tmax = 250
env.change = data.frame(time=seq(1,tmax,by=0.1),env.dynamic=9)

set.seed(4) # pick some useful trait values
trait_values = sort(runif(n=nsp))
trait_values[nsp] = 1.0 # last species gets biggest trait value
t1 = 25
t2 = t1+1 # not used

set.seed(1)
N_t0 = c(rep(0.01,nsp-1),0) #all species but one
env_change_response_ptA <- predict_environmental_change_response(nsp = nsp,
                                                                 env.change = env.change %>% filter(time <= t1),
                                                                 trait.distribution = 'forced',
                                                                 trait_values = trait_values,
                                                                 mechanism = "niche difference non normalized",
                                                                 Nmin = 0.001,
                                                                 initial.abundance = N_t0,
                                                                 growth.rate = 0.5,
                                                                 extinction = T,
                                                                 plot.species.dynamics = F,
                                                                 plot.response.diagram = F)

n_t1 = env_change_response_ptA$comm_dynamic %>% filter(time==t1) %>% pull(n)

set.seed(1)
env_change_response_ptB <- predict_environmental_change_response(nsp = nsp,
                                                                 env.change = env.change %>% filter(time > t1 & time <= t2),
                                                                 trait.distribution = 'forced',
                                                                 trait_values = trait_values,
                                                                 mechanism = "niche difference non normalized",
                                                                 Nmin = 0.001,
                                                                 initial.abundance = n_t1,
                                                                 growth.rate = 0.5,
                                                                 extinction = T,
                                                                 plot.species.dynamics = F,
                                                                 plot.response.diagram = F)

n_t2 = env_change_response_ptB$comm_dynamic %>% filter(time==t2) %>% pull(n)

set.seed(1)
env_change_response_ptC <- predict_environmental_change_response(nsp = nsp,
                                                                 env.change = env.change %>% filter(time > t2),
                                                                 trait.distribution = 'forced',
                                                                 trait_values = trait_values,
                                                                 mechanism = "niche difference non normalized",
                                                                 Nmin = 0.001,
                                                                 initial.abundance = n_t2,
                                                                 growth.rate = 0.5,
                                                                 extinction = T,
                                                                 plot.species.dynamics = F,
                                                                 plot.response.diagram = F)

env_change_all = rbind(env_change_response_ptA$comm_dynamic, env_change_response_ptB$comm_dynamic, env_change_response_ptC$comm_dynamic)







set.seed(1)
n_t1_perturbed = n_t1
n_t1_perturbed[nsp] = 0.01 # add species #10 at low density
env_change_response_ptB_perturbed <- predict_environmental_change_response(nsp = nsp,
                                                                           env.change = env.change %>% filter(time > t1 & time <= t2),
                                                                           trait.distribution = 'forced',
                                                                           trait_values = trait_values,
                                                                           mechanism = "niche difference non normalized",
                                                                           Nmin = 0.001,
                                                                           initial.abundance = n_t1_perturbed,
                                                                           growth.rate = 0.5,
                                                                           extinction = T,
                                                                           plot.species.dynamics = F,
                                                                           plot.response.diagram = F)

n_t2_perturbed = env_change_response_ptB_perturbed$comm_dynamic %>% filter(time==t2) %>% pull(n)

set.seed(1)
env_change_response_ptC_perturbed <- predict_environmental_change_response(nsp = nsp,
                                                                           env.change = env.change %>% filter(time > t2),
                                                                           trait.distribution = 'forced',
                                                                           trait_values = trait_values,
                                                                           mechanism = "niche difference non normalized",
                                                                           Nmin = 0.001,
                                                                           initial.abundance = n_t2_perturbed,
                                                                           growth.rate = 0.5,
                                                                           extinction = T,
                                                                           plot.species.dynamics = F,
                                                                           plot.response.diagram = F)

env_change_all_perturbed = rbind(env_change_response_ptA$comm_dynamic, env_change_response_ptB_perturbed$comm_dynamic, env_change_response_ptC_perturbed$comm_dynamic)




# prep data for plotting
eca.trait.unperturbed = env_change_all %>% 
  mutate(label=LETTERS[factor(trait.i,levels=sort(unique(trait.i)))])
eca.trait.perturbed = env_change_all_perturbed %>% 
  mutate(label=LETTERS[factor(trait.i,levels=sort(unique(trait.i)))])

g_addition_unperturbed = ggplot(eca.trait.unperturbed, aes(y = n, x = time)) +
  geom_line(aes(col = trait.i, group = i),size=1) +
  scale_color_viridis("Trait value",limits=c(0,1)) +
  theme_bw()+
  theme(legend.position = "none") +
  #scale_y_sqrt() +
  xlab("Time (t)") + ylab("Abundance (N)") +
  geom_label_repel(data=eca.trait.unperturbed %>% filter(time >= (tmax-0.1) & time < tmax),aes(x=time,y=n,label=label,color=trait.i),max.overlaps=Inf) +
  #geom_label_repel(data=env_change_all %>% filter(time==tmax),aes(x=time,y=n,label=label,color=trait.i)) +
  ggtitle('Uncontrolled dynamics')


g_addition_perturbed = ggplot(eca.trait.perturbed, aes(y = n, x = time)) +
  geom_line(aes(col = trait.i, group = i),size=1) +
  scale_color_viridis("Trait value",limits=c(0,1)) +
  theme_bw()+
  theme(legend.position = "none") +
  #scale_y_sqrt() +
  xlab("Time (t)") + ylab("Abundance (N)") +
  geom_label_repel(data=eca.trait.perturbed %>% filter(time==tmax),aes(x=time,y=n,label=label,color=trait.i),max.overlaps=Inf) +
  #geom_label_repel(data=env_change_all %>% filter(time==tmax),aes(x=time,y=n,label=label,color=trait.i)) +
  ggtitle('Controlled dynamics')


introductions = data.frame(time=seq(0,tmax,by=0.1),deltaN=0)
introductions$deltaN[introductions$time==t1] = 0.01

g_addition_controlling_variable = ggplot(introductions,aes(x=time,y=deltaN)) +
  geom_line(size=1) +
  geom_hline(yintercept = 0,size=1,linetype='dotted',color='gray') +
  theme_bw() +
  xlab("Time (t)") + ylab("∆N of species J") +
  ggtitle('Controlled variable')


# look at trait values
eca.trait.unperturbed %>% select(i, label, trait.i) %>% unique

g_addition = ggarrange(g_addition_controlling_variable, g_addition_unperturbed, g_addition_perturbed,
                       align='hv',common.legend = TRUE, legend='bottom',labels='auto',nrow=1,ncol=3)
ggsave(g_addition, file='g_control_addition.png',width=9,height=4)











