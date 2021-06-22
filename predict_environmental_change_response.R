
################################################################################
################### Predict environmental change response ######################
########################################### Gauzere P. and Benjamin Blonder#####

# The function predict the temporal dynamic of a given community under temporal environmental change
# simulated with function simulate_environmental_change()

# the difference with previous functions is that the A matrix is recalculated at each time step based on the environmental value at the time t

# function argument are  : 
##same as predict_demographic_model_parameters :
  # nsp : number of species in the community before assembly [numeric]
  #  trait.distribution (for now only one possible)
  #  mechanism (for now only one possible)

##same as predict_demographic_model_parameters : 
  #  Nmin
  #  initial.abundance 
  #  growth.rate (note that the invader will have r = trait value)
  #  extinction 

##new arguments 
  # env.change [data.frame] : the data frame returned by simulate_environmental_change()
  # plot.species.dynamics [logical] : if TRUE, plot the community dynamic along with environment dynamic
  # traits_to_test [numeric] : if TRUE, plot the community weighted mean of the trait, the with environment dynamic, an the community response diagram  


### load library
require(tidyverse)
require(deSolve)
require(rootSolve)
require(viridis)
require(ggpubr)

calculate_jacobian_glv <- function(r_vec, A_mat, N_vec, allow_extinction=TRUE, threshold.abundance=0.005)
{
  #J_mat = matrix(NA,nnrow=nrow(A_mat),ncol=ncol(A_mat))
  
  N_mat = matrix(data=N_vec,nrow=nrow(A_mat),ncol=ncol(A_mat))
  
  #print(N_mat)
  
  r_mat = diag(r_vec)
  
  #print(r_mat)
  
  J_mat = r_mat + A_mat %*% N_mat
  
  if (allow_extinction==TRUE)
  {
    # jacobian is 0 with respect to any variable if a species is extinct
    which_species = which(N_vec < threshold.abundance)
    J_mat[which_species,] <- 0
  }
  
  #print(J_mat)
  
  return(J_mat)
}

predict_environmental_change_response <- function(nsp = 10,
                                                  env.change,
                                                  trait.distribution = "uniform",
                                                  trait_values = NULL,
                                                  mechanism = "niche difference",
                                                  Nmin = 0.001,
                                                  initial.abundance = 0.01,
                                                  growth.rate = 0.2,
                                                  extinction = F,
                                                  plot.species.dynamics = T,
                                                  plot.response.diagram = T, 
                                                  probs.niche.difference=0.1) {
  
  
  
  # we take the timeseries from env.change
  time <- env.change$time
  
  
  source("predict_demographic_model_parameters.R")
  alpha.df <- predict_demographic_model_parameters(rep = 1,
                                                   nsp = nsp,
                                                   trait.distribution = trait.distribution,
                                                   trait_values = trait_values,
                                                   mechanism = mechanism, 
                                                   env = env.change$env.dynamic,
                                                   probs.niche.difference = probs.niche.difference
                                                     )
  # we add an epxlicit time columns to alpha.df
  alpha.df$time <-
    rep(time,
        each = nsp ^ 2 * length(mechanism) * length(trait.distribution))
  
# classic now... define No...
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
  
  
  
  # ... and r given the arguments
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
      r = unique(alpha.df$trait.i)
    }
  }
  
  #define equation system
  #note that this is different than before. For each t we have a new A matrix, which depend on environment value at time t
  eqs <- function(t, x, params) {
    A <-
      xtabs(interaction.coef ~ i + j, data = alpha.df[alpha.df$time == ceiling(t), ])
    x[x < Nmin] = 0
    dx = x * (r + A %*% x) # add a constant here
    
    # dx[x<Nmin]=0 # extinct species don't grow anymore - remove to allow for reinvasions
    return(list(dx))
  }
  

  print("numerical simulation running...")
  result = ode(y = No, times = time, func = eqs)
  print("...done")
  
  # calculate the jacobian matrix
  empty_list = vector(mode="list",length=length(time))
  r_list = empty_list
  A_list = empty_list
  N_list = empty_list
  J_list = empty_list
  
  for (i in 1:length(time))
  {
    # copy all the r values (no temporal change)
    r_list[[i]] = r
    # update the A values
    A_list[[i]] = xtabs(interaction.coef ~ i + j, data = alpha.df[alpha.df$time == ceiling(time[i]), ])
    N_list[[i]] = as.numeric(result[i,2:ncol(result)]) # copy out the abundances - dropping 1st column (time)
    J_list[[i]] = calculate_jacobian_glv(r_list[[i]], A_list[[i]], N_list[[i]])
  }
  
  detJ_timeseries = data.frame(time=time,detJ=sapply(J_list,det))

  #built comm_dynamic dataframe
  comm_dynamic <-
    gather(as.data.frame(result),
           key = "i",
           value = "n",
           -time)
  comm_dynamic$i <- as.integer(comm_dynamic$i)
  comm_dynamic <- left_join(comm_dynamic, unique(alpha.df[, c("i", "trait.i")]))
  

  #compute community indices at each time
  cwm_dynamic <-
    left_join(comm_dynamic, unique(alpha.df[, c("i", "trait.i")])) %>%
    group_by(time) %>%
    summarize("n_tot" = sum(n),
              "cwm" = mean( trait.i   * (n/n_tot)),
               "cwv" = sum ( trait.i^2 * (n/n_tot)) - cwm) 
  cwm_dynamic$env <- env.change$env.dynamic
  
  
  #plot species dynamics, if needed
  if (plot.species.dynamics == T) {
    require(cowplot)
    env.plot <-  
      ggplot(cwm_dynamic) +
      geom_line(aes(x = time, y = env), col = "brown4", size = 1) +
      scale_color_viridis() +
      theme_classic()
    
    comm.plot <- 
      ggplot(comm_dynamic, aes(y = n, x = time)) +
      geom_line(aes(col = trait.i, group = i), size = 1) +
      scale_color_viridis() +
      theme_classic()+
      theme(legend.position = "bottom")
    
    
    print(plot_grid(env.plot, comm.plot, ncol = 1, align = "v", rel_heights = c(1,2)))
  }
  
  #plot community dynamics and response diagram, if needed
  if (plot.response.diagram == T) {
    require(cowplot)
    cwm.plot <-  
      ggplot(cwm_dynamic) +
      geom_line(aes(x = time, y = cwm), col = "darkgreen", size = 1) +
      scale_color_viridis_d() +
      theme_classic()
    
    response.diagram <-
      ggplot(cwm_dynamic) +
      geom_path(aes(x = env, y = cwm, col = time), size = 1) +
      scale_color_viridis() +
      theme_classic()+
      theme(legend.position = "bottom")
    
    
    print(plot_grid(env.plot, cwm.plot, response.diagram, ncol = 1, rel_heights = c(1, 1, 3), align = "v"))
  }
#return both dataframes
cwm = left_join(comm_dynamic, cwm_dynamic)
return(list("comm_dynamic"=cwm, detJ=detJ_timeseries, J_list=J_list))
}
