################################################################################
##### Predict demographic model parameter (from trait x environment) ###########
########################################### Gauzere P. #########################

#note : you can use this script as a function, 
#or you can uncomment the next row (#33) and the last row (#162) 
#to run the script step by step for better understanding and tuning possibilities 
#and use the section #1 define parameters

# requires tidyverse to be installed 

# function argument are : 
  # rep : number of replicats [numeric] 
  # nsp : number of species in the community before assembly [numeric]
  # env : environmental/stress gradient. Only positive values [0, +Inf] [numeric]
  # trait.distribution : distribution of functional traits [character string]
                  #choices are  : "gaussian", "uniform", "poisson", "bimodal"
                  #possibility to put any combinations. e.g : c("gaussian", "uniform"), c("gaussian", "uniform", "poisson", "bimodal")
  # mechanism : assembly mechanism hypothesis [character string]
                  #choices are : "niche difference", "competitive dominance","extreme facilitation"
                  # possibility to put any combinations. 
                  

# return a dataframe with columns : 
  # i : index of species i                 
  # j : index of species j
  # env : environmental/stress value/gradient/forcing
  # mechanism : assembly mechanism used
  # trait.distribution : trait distribution used
  # replicat : replicate index
  # interaction.coef : the value of interaction coefficient
  # trait.i : trait value of species i          
  # trait.j : trait value of species j

predict_demographic_model_parameters <- function(rep = 10,
                                                 nsp = 10,
                                                 env = seq(from = 0, to = 10, by = 0.1),
                                                 trait.distribution = "uniform",
                                                 mechanism = "competitive dominance", 
                                                 trait_values = seq(0, 1, length.out = nsp)) {
  
  require(tidyverse)

#### 1 load funtional form under assembly hypotheses 
source("compute_interaction_from_niche_difference.R")
source("compute_interaction_from_competitive_dominance.R")
source("compute_interaction_from_extreme_facilitation.R")# this one often brings unstability when modelizing the community dynamics
  
##### 3 prepare empty matrices and dataframe given parameters ####
# create empty matrices number and size will depend on parameters. 
alpha <-
  array(
    rep(NA, nsp ^ 2),
    dim = c(
      nsp,
      nsp,
      length(env),
      length(mechanism),
      length(trait.distribution),
      rep
    ),
    dimnames = list(1:nsp, 1:nsp, env, mechanism, trait.distribution, 1:rep)
  )

# reshape matrix to long dataframe (for computation purpose)
alpha.df <- reshape2::melt(alpha)
colnames(alpha.df) <-
  c("i",
    "j",
    "env",
    "mechanism",
    "trait.distribution",
    "replicat",
    "interaction.coef")
alpha.df$trait.i <- NA
alpha.df$trait.j <- NA

#### 4 compute interaction coefficient given parameters for each replicat####

for(i in 1:rep){ #looping trhough replicats
sub_alpha.df <- alpha.df[alpha.df$replicat == i,]  

if(trait.distribution == "forced"){
  sub_alpha.df$trait.i <- rep(trait_values, length.out = nrow(sub_alpha.df))
  sub_alpha.df$trait.j <- rep(trait_values, each = nsp)
} else{

## random generation of traits values following distributions
#gaussian
traits_norm    <- rnorm(n = nsp)
#uniform
traits_uniform <- runif(n = nsp)
#poisson
traits_poisson <- rpois(n = nsp, lambda = 1)
#bimodal (normal mixture)
y1 = rnorm(nsp, -1, 0.5)
y2 = rnorm(nsp, 1, 0.5)
w  = rbinom(nsp, 1, .5) # 50:50 random choice
traits_bimod <- as.numeric(scale(w*y1 + (1-w)*y2)) # normal mixture

## add i and j trait value depending on distributions in sub_alpha.df
# i
sub_alpha.df$trait.i[sub_alpha.df$trait.distribution == "gaussian"] <- traits_norm
sub_alpha.df$trait.i[sub_alpha.df$trait.distribution == "uniform"]  <- traits_uniform
sub_alpha.df$trait.i[sub_alpha.df$trait.distribution == "poisson"]  <- traits_poisson
sub_alpha.df$trait.i[sub_alpha.df$trait.distribution == "bimodal"]  <- traits_bimod

#j
sub_alpha.df$trait.j[sub_alpha.df$trait.distribution == "gaussian"] <- rep(traits_norm, each=nsp)
sub_alpha.df$trait.j[sub_alpha.df$trait.distribution == "uniform"]  <- rep(traits_uniform, each=nsp)
sub_alpha.df$trait.j[sub_alpha.df$trait.distribution == "poisson"]  <- rep(traits_poisson, each=nsp)
sub_alpha.df$trait.j[sub_alpha.df$trait.distribution == "bimodal"]  <- rep(traits_bimod, each=nsp)
}

# compute interaction coeffcients depending on mechanisms and environment value
sub_alpha.df$interaction.coef[sub_alpha.df$mechanism == 'niche difference'] <- 
  niche_difference(x = sub_alpha.df[sub_alpha.df$mechanism == 'niche difference',])

sub_alpha.df$interaction.coef[sub_alpha.df$mechanism == 'competitive dominance'] <- 
  competitive_dominance(x = sub_alpha.df[sub_alpha.df$mechanism == 'competitive dominance',] )

sub_alpha.df$interaction.coef[sub_alpha.df$mechanism == 'extreme facilitation'] <- 
  extreme_facilitation( x= sub_alpha.df[sub_alpha.df$mechanism == 'extreme facilitation',])

alpha.df[alpha.df$replicat == i,] <- sub_alpha.df
rm(sub_alpha.df)

print(paste("replicat", i,"on", rep,"| env", env, sep = " "))  
}
return(alpha.df)}
