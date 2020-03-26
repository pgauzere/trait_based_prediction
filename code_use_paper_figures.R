library(viridis)

source("predict_demographic_model_parameters.R")

## we simulate a continuous environmental gradient under the competitive dominance
 #read the forewords of the function predict_demographic_model_parameters()
alpha.df <- predict_demographic_model_parameters(
  rep = 10,
  nsp = 10,
  env = seq(from = 0, to = 10, by = 0.1),
  trait.distribution = "uniform",
  mechanism = "competitive dominance")

##### Figure 4 #####
### a interaction coeff distribution 
ggplot(alpha.df%>% filter( i != j), #remove intraspecific coefficient for plot
       aes(x=interaction.coef, fill=env, group = env))+
  geom_histogram()+
  theme_classic()+
  theme(legend.position="bottom")+
  scale_fill_viridis()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), trans= "sqrt")

### we can play with different mechanisms and distributions
# predict coefficient under multiple hypotheses and trait distributions 
alpha.df <- predict_demographic_model_parameters(
  rep = 10,
  nsp = 10,
  env = seq(from = 0, to = 10, by = 0.1),
  trait.distribution = c("uniform","gaussian", "bimodal"),
  mechanism = c("niche difference", "competitive dominance", "extreme facilitation")
)

#plot all of it to compare
alpha.df %>% filter( i != j) %>% #remove intraspecific coefficient
  ggplot(aes(x=interaction.coef, fill=env, group = env))+
  geom_histogram()+
  theme_classic()+
  theme(legend.position="bottom")+
  scale_fill_viridis()+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0), trans= "sqrt")+
  facet_wrap(~ mechanism + trait.distribution)


### b-c plot netwrok #####
source("predict_coexistence_outcome.R")

# we can use predict_demographic_model_parameters() with only 2 environmental values
# choose environmental value : env = 0 for fig 4.b, env = 10 for fig 4.c
alpha.df <- predict_demographic_model_parameters(
  rep = 1,
  nsp = 10,
  env = c(0, 10),
  trait.distribution = "uniform",
  mechanism = "competitive dominance")

#then we use the function predict_coexistence_outcome() to plot the netwrok with "plot.netwrok=T"
#please read the forewords of the function predict_coexistence_outcome() 
alpha.df %>% 
  filter(env == 0) %>%
  predict_coexistence_outcome(plot.network = T)

alpha.df %>% 
  filter(env == 10) %>%
  predict_coexistence_outcome(plot.network = T)

# note that color scale in paper Figure 4 is adjusted limits[-1,0] to facilitate comparison between network 

##### Figure 5 #####

# we need a good environmental definition and trait coverage to interpolate coexistence outcome surface
alpha.df <- predict_demographic_model_parameters(
  rep = 100,
  nsp = 10,
  env = seq(from = 0, to = 10, by = 0.1))

# WARNING computation can take a bit of time. To avoid that : 
# load("alpha.properties_10sp_100_replicates_100env.Rdata")

alpha.properties <-
    alpha.df %>%
    group_by(env, replicat) %>%
    group_modify( ~ predict_coexistence_outcome(.x))

save(alpha.properties, file = "alpha.properties_10sp_100_replicates_100envBIS.Rdata")


ggplot(alpha.properties, aes(x = env, y = mean_trait))+
  geom_point(aes(col = feasibility))+
  scale_color_viridis_d()

ggplot(alpha.properties, aes(x = env, y = mean_trait))+
  geom_point(aes(col = stability))+
  scale_color_viridis_d()


ggplot(alpha.properties, aes(x = env, y = var_trait))+
  geom_point(aes(col = n_realized))+
  scale_color_viridis()

ggplot(alpha.properties, aes(x = env, y = var_trait))+
  geom_point(aes(col = cwv - var_trait))+
  scale_color_viridis()

ggplot(alpha.properties, aes(x = env, y = mean_trait))+
  geom_point(aes(col = cwm - mean_trait))+
  scale_color_viridis()

## feasability probability
library(akima)
feasibility_interpolation <-
  interp(
    alpha.properties$env,
    alpha.properties$var_trait,
    alpha.properties$feasibility,
    nx = 100,
    ny = 100,
    duplicate = "strip"
  )
li.zmin <- min(feasibility_interpolation$z,na.rm=TRUE)
li.zmax <- max(feasibility_interpolation$z,na.rm=TRUE)
breaks <- pretty(c(li.zmin,li.zmax),10)
colors <- viridis(length(breaks)-1)
with(feasibility_interpolation, image (x,y,z, breaks=breaks, col=colors, ylim = c(0.05, 0.12)))

##stability
stab_interpolation <- interp(alpha.properties$env, alpha.properties$mean_trait, alpha.properties$stability, nx=100, ny=100, duplicate = "strip")
li.zmin <- min(stab_interpolation$z,na.rm=TRUE)
li.zmax <- max(stab_interpolation$z,na.rm=TRUE)
breaks <- pretty(c(li.zmin,li.zmax),10)
colors <- viridis(length(breaks)-1)
with(stab_interpolation, image (x,y,z, breaks=breaks, col=colors, ylim = c(0.3, 0.7)))
# contour(stab_interpolation, levels=breaks, add=TRUE)


##n_realized
n_interpolation <- interp(alpha.properties$env, alpha.properties$mean_trait, alpha.properties$n_realized, nx=100, ny=100, duplicate = "strip")
li.zmin <- min(n_interpolation$z,na.rm=TRUE)
li.zmax <- max(n_interpolation$z,na.rm=TRUE)
breaks <- pretty(c(li.zmin,li.zmax),10)
colors <- viridis(length(breaks)-1)
with(n_interpolation, image (x,y,z, breaks=breaks, col=colors, ylim = c(0.3, 0.7)))
 contour(n_interpolation, levels=breaks, add=TRUE)


##### Figure  6 #####

##### Figure  7 #####
# here we are interested of the temporal dynamic of a given community
# in this case the environmental gradient is explicitely temporal
# we use the function predict_environmental_change_response(), which inherits from
# predict_model_parameters() and predict_coexistence_outcome()
# please read the forewords of the function predict_coexistence_outcome() 



