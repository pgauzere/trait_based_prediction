
compute_invasibility_surface <- function(
  nsp = 10,
  env = seq(from = 0, to = 10,  by = 0.5),  
  trait.distribution = "uniform",  
  mechanism = "competitive dominance" ,
  Nmin = 0.0001 ,
  initial.abundance = 0.0005 ,
  growth.rate = 0.5,
  n.time.step = 100 ,
  extinction = F ,
  plot.invasion = F,
  traits_to_test = seq(0,1,0.05)){
  
  # 
  # nsp = 4
  # env = seq(from = 0, to = 10,  by = 1)
  # trait.distribution = "uniform"
  # mechanism = "competitive dominance"
  # Nmin = 0.0001
  # initial.abundance = 0.0005
  # growth.rate = 0.5
  # n.time.step = 100
  # extinction = F
  # plot.invasion = F
  # traits_to_test = seq(0,1,0.1)
  # 

  
  ### load library
  require(tidyverse)
  require(deSolve)
  require(rootSolve)
  require(igraph)
  require(ggraph)
  require(NetIndices)
  require(viridis)

  source("create_empty_alpha.df.R")
  source("predict_demographic_model_parameters.R")
  source("predict_coexistence_outcome.R")
  source("compute_interaction_from_competitive_dominance.R")
  source("compute_interaction_from_niche_difference.R")
  source("compute_interaction_from_extreme_facilitation.R")
  
  
alpha.df <- predict_demographic_model_parameters(rep = 1, 
                                  nsp = nsp, 
                                  env = env, 
                                  trait.distribution = trait.distribution, 
                                  mechanism = mechanism)

# build a new empty alpha.df with one more species (the invader)
alpha.df.inv <- create_empty_alpha.df(
                      rep = 1, 
                      nsp = nsp+1, 
                      env = env, 
                      trait.distribution =trait.distribution, 
                      mechanism = mechanism)

# fill alpha.df.inv with known interaction coefficients from alpha.df. all but invader species 
alpha.df.inv[alpha.df.inv$i != nsp+1 & alpha.df.inv$j != nsp+1, c("interaction.coef", "trait.i", "trait.j")] <- alpha.df[, c("interaction.coef", "trait.i", "trait.j")]
alpha.df.inv[alpha.df.inv$i != nsp+1 & alpha.df.inv$j == nsp+1, "trait.i"] <- rep(unique(alpha.df$trait.i))
alpha.df.inv[alpha.df.inv$i == nsp+1 & alpha.df.inv$j != nsp+1, "trait.j"] <- rep(unique(alpha.df$trait.j))

#create empty matrix for invasibility surface
res_invasion <-
  matrix(
    data = NA,
    nrow = length(env),
    ncol = length(traits_to_test),
    dimnames = list(as.character(env),
                    as.character(traits_to_test))
  )

#### 1 simulate the community dynamic for a given environment, without invader ####
# the idea is to reach equilibrium under the environemtal conditions

for(i in 1:length(env)){
  
  x <- alpha.df[alpha.df$env == env[i],]

  ### define ODE equations system for multi-species dynamics
  # initial abundance
  if (class(initial.abundance) == "numeric") {
    if (length(initial.abundance) == 1) {
      No = rep(initial.abundance, nsp)
    } else{
      No = initial.abundance
    }
  } else {
    if (initial.abundance == "random") {
      No = runif(nsp)
    }
  }
  
  # growth rate
  if (class(growth.rate) == "numeric") {
    if (length(growth.rate) == 1) {
      r = rep(growth.rate, nsp)
    } else{
      r = growth.rate
    }
  } else {
    if (growth.rate == "random") {
      # intrinsic growth rates sampled from normal distribution
      r = rnorm(nsp, 1, 0.1)
    }
    if (growth.rate == "trait") {
      r = x$trait.i[1:nsp]
    }
  }
  
  ##define equation system
  eqs <- function(t, n, params) {
    n[n < Nmin] = 0
    dn = n * (r + A %*% n) # add a constant here
    if(extinction == T){dn[n<Nmin]=0} else{} # extinct species don't grow anymore - remove to allow for reinvasions
    return(list(dn))
  }  
  
  #convert table into matrix
  A <- xtabs(interaction.coef ~ i + j, data=x)
  
  ### Numerical simulations

  #run ODE dynamics
  result = ode(y = No, times=1:n.time.step, func=eqs)
  comm_final <- as.numeric(result[n.time.step,-1])
  
  ##### 2 we invade the community sequentially with each trait value #### 
  
  for(j in 1:length(traits_to_test)){
    
    # we compute the interaction coefficients of invders wiht resident community under environment
    alpha.df.inv[alpha.df.inv$i == nsp+1, "trait.i"] <- traits_to_test[j]
    alpha.df.inv[alpha.df.inv$j == nsp+1, "trait.j"] <- traits_to_test[j]
    
    # if(mechanism == "competitive_dominance") {
      alpha.df.inv[alpha.df.inv$i == nsp + 1 |
                     alpha.df.inv$j == nsp + 1, "interaction.coef"] <-
        competitive_dominance(x = alpha.df.inv[alpha.df.inv$i == nsp + 1 |
                                                 alpha.df.inv$j == nsp + 1, ])
    # } else{
    # if (mechanism == "niche difference") {
    #   alpha.df.inv[alpha.df.inv$i == nsp + 1 |
    #                  alpha.df.inv$j == nsp + 1, "interaction.coef"] <-
    #     niche_difference(x = alpha.df.inv[alpha.df.inv$i == nsp + 1 |
    #                                         alpha.df.inv$j == nsp + 1, ])
    # } else{
    #   alpha.df.inv[alpha.df.inv$i == nsp + 1 |
    #                  alpha.df.inv$j == nsp + 1, "interaction.coef"] <-
    #     extreme_facilitation(x = alpha.df.inv[alpha.df.inv$i == nsp + 1 |
    #                                             alpha.df.inv$j == nsp + 1, ])
    # }}
    # 
    x2 <- alpha.df.inv %>% filter(env == env[i])
    
    #adding initial abundance of invader. Invader is introduced at low density, i.e Nmin *5
    No2 = c(comm_final, Nmin*100) # introduction at low density

    #convert table into matrix
    A2 <- xtabs(interaction.coef ~ i + j, data=x2)
    
       #  adding growth rate invader
    # if (class(growth.rate) == "numeric") {
    #   if (length(growth.rate) == 1) {
    #     r2 = c(r[1:nsp], growth.rate)
    #   } else{
    #     r2 = c(r, growth.rate)
    #   }
    # } else {
    #   if (growth.rate == "random") {
    #     # intrinsic growth rates sampled from normal distribution
    #     r2 = c(r[1:nsp], rnorm(1, 1, 0.1))
    #   }
      # if (growth.rate == "trait") {
        r2 = c(r[1:nsp], traits_to_test[j])
    #   }
    # }
    
    ##define equation system
    eqsInv <- function(t, n, params) {
      # n[n < Nmin] = 0
      dn = n * (r2 + A2 %*% n) # add a constant here
      if(extinction == T){dn[n<Nmin]=0} else{} # entinct species don't grow anymore - remove to allow for reinvasions
      return(list(dn))
    } 
    
    #run numeric simulations with resident + invader
    result = ode(y = No2, times=1:n.time.step, func=eqsInv)
    
    
    # plot invasion dynamic
    comm_dynamic <- gather(as.data.frame(result), key="species_i", value="n", -time)
    if(plot.invasion == T){
    print(
      ggplot(comm_dynamic, aes(y=n, x=time, col=species_i))+
        geom_line()+
        geom_line(data = comm_dynamic[comm_dynamic$species_i == as.character(nsp+1),], aes(y=n, x=time), col="red")+
        scale_color_viridis_d(guide=F)+
        geom_hline(yintercept = 0.0001, linetype = 2)+
        theme_classic()+
        labs(title = paste("trait value =",traits_to_test[j]))+
        ylim(0,1)
    )
    }
    #store results of invasion analysis 
    res_invasion[i,j] <- comm_dynamic[nrow(comm_dynamic), "n"]
    
    #monitor progress
    print(paste("invasion : ","env[", i,"]/",length(unique(alpha.df$env)), "; trait[", j, "]/", length(traits_to_test),  "-- done", sep =''))
  }
}
return(res_invasion)
}

