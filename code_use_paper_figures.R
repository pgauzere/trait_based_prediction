################################################################################
##### Use the code and re-produce figures from the Article######################
##### Trait-based prediction and control: a roadmap for community ecology  #####
########################################### Gauzere P. #########################

library(viridis)

############### A Predicting demographic model parameters ######################

#We will use the function predict_model_parameters()
#this basically implements the idea presented in chapter "Predicting demographic model parameters"
# and Box 2b
#load the function
source("predict_demographic_model_parameters.R")
#read the forewords of the function to know how to use it

## we simulate a continuous environmental gradient under the competitive dominance
alpha.df <- predict_demographic_model_parameters(
  rep = 10,
  nsp = 10,
  env = seq(from = 0, to = 10, by = 0.1),
  trait.distribution = "uniform",
  mechanism = "competitive dominance")

##### Figure 4a #####
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


############### B Scaling up to coexistence outcomes ######################

# We will use the function predict_coexistence_outcome()
# this basically implements the idea presented in chapter "Scaling up to coexistence outcomes"
# and Box 2c
# load the function
source("predict_coexistence_outcome.R")
# read the forewords of the function to know how to use it


### Figure 4b-c : plot netwrok #####

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


#### the function compute a lot of network features and coexistence, but is time demanding 
# note that the resolution of the ODE sometimes won't be possible if you choose weird parameters
# keep it easy
# let's play
# we need a good environmental definition and trait coverage to interpolate coexistence outcome surface
alpha.df <- predict_demographic_model_parameters(
  rep = 20,
  nsp = 6,
  env = seq(from = 0, to = 10, by = 0.2))


# for computation purpose, we use tidyR and group_by() - group_modify(), but 
# you can use the function to loop on the different subset of alpha.df if you prefer

alpha.properties <-
    alpha.df %>%
    group_by(env, replicat) %>%
    group_modify( ~ predict_coexistence_outcome(.x,
                                                Nmin = 0.0001,
                                                initial.abundance = 0.01,
                                                growth.rate = 0.5,
                                                n.time.step = 500,
                                                extinction = T,
                                                plot.dynamic = T,
                                                plot.network = T, 
                                                network.threshold = 0.05
                                                ))

ggplot(alpha.properties, aes(x = env, y = mean_trait))+
  geom_point(aes(col = feasibility, size = n_realized, alpha = n_realized))+
  scale_color_viridis_d()+
  theme_classic()

ggplot(alpha.properties, aes(x = env, y = mean_trait))+
  geom_point(aes(col = stability, size = n_realized, alpha = n_realized))+
  scale_color_viridis_d()+
  theme_classic()

ggplot(alpha.properties, aes(x = env, y = var_trait))+
  geom_point(aes(col = n_realized, size = n_realized, alpha = n_realized))+
  scale_color_viridis()+
  theme_classic()

ggplot(alpha.properties, aes(x = env, y = var_trait))+
  geom_point(aes(col = cwv , size = n_realized, alpha = n_realized))+
  scale_color_viridis()+
  theme_classic()

ggplot(alpha.properties, aes(x = env, y = mean_trait))+
  geom_point(aes(col = cwm - mean_trait))+
  scale_color_viridis()

##### Figure 5 #####
# first , we need a lot of replicates to sample a variation in trait mean due to randomness
alpha.df <- predict_demographic_model_parameters(
  rep = 100,
  nsp = 10,
  env = seq(from = 0, to = 10, by = 0.1))

#then we can compute the coexistence outcome, but it can take a lot of time
# to avoid that, you can load the output directly :
load("alpha.properties_10sp_100_replicates_100env.Rdata")

# alpha.properties <-
#   alpha.df %>%
#   group_by(env, replicat) %>%
#   group_modify( ~ predict_coexistence_outcome(.x))


# To create figure 5, we interpolate the predicted over the full surface of trait - env

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

############### C 
 # We will use the function predict_invasibility()
 # this basically implements the idea presented in chapter "Predicting invasibility from trait×environment interactions"
 # load the function
 source("predicting_invasibility.R")
 # read the forewords of the function to know how to use it
res_invasion <-
   predict_invasibility_surface(
     nsp = 10,
     env = seq(from = 0, to = 10,  by = 1),
     trait.distribution = c("uniform"),
     mechanism = "competitive dominance",
     Nmin = 0.0001 ,
     initial.abundance = 0.0005 ,
     growth.rate = 0.5,
     n.time.step = 250,
     extinction = T,
     plot.invasion = T,
     traits_to_test = seq(0,1, 0.1)
   )

ggplot(reshape2::melt(res_invasion) 
        , aes(Var1, Var2, fill = value))+
   geom_tile()+
   scale_fill_viridis(limits = c(0,1))+
   theme_minimal()+
   labs(x = "environment", y = "trait value", fill = "invasibility")+
   scale_x_continuous(expand=c(0,0)) +
   scale_y_continuous(expand=c(0,0))
 
##### Figure  6 #####

 
############### D Predicting environmental change response

# here we are interested of the temporal dynamic of a given community
# in this case the environmental gradient is explicitely temporal
# we use the function predict_environmental_change_response(), which inherits from
# predict_model_parameters() and predict_coexistence_outcome()
# please read the forewords of the function predict_coexistence_outcome() 

 ##### Figure  7 #####

# We will use the function predict_invasibility()
# this basically implements the idea presented in chapter "Predicting invasibility from trait×environment interactions"
# load the function
source("predict_environmental_change_response.R")
# read the forewords of the function to know how to use it

test <- predict_environmental_change_response()
test

