


predict_environmental_change_response <- function(nsp = 10,
                                                  env.change = "sinusoidal",
                                                  trait.distribution = "uniform",
                                                  mechanism = "niche difference",
                                                  Nmin = 0.001,
                                                  initial.abundance = 0.01,
                                                  growth.rate = 0.2,
                                                  n.time.step = 100,
                                                  extinction = F,
                                                  plot.species.dynamics = T,
                                                  plot.response.diagram = T) {
  
  # nsp = 20
  # env.change = "sinusoidal"
  # trait.distribution = "uniform"
  # mechanism = "niche difference"
  # Nmin = 0.001
  # initial.abundance = 0.01
  # growth.rate = 0.2
  # n.time.step = 100
  # extinction = F
  # plot.species.dynamics = T
  # plot.response.diagram = T
  
  time <- 1:n.time.step
  
  source("predict_demographic_model_parameters.R")
  alpha.df <- predict_demographic_model_parameters(rep = 1,
                                                   nsp = nsp,
                                                   trait.distribution = trait.distribution,
                                                   mechanism = mechanism, 
                                                   env =
                                                     if (env.change == "sinusoidal") {
                                                       1 + sin(time / (n.time.step / 10))
                                                     } else{
                                                       time / 10
                                                     })
  
  alpha.df$time <-
    rep(time,
        each = nsp ^ 2 * length(mechanism) * length(trait.distribution))
  
  # plot(alpha.df$time, alpha.df$env)
  
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
      r = unique(alpha.df$trait.i)
    }
  }
  
  #define equation system
  eqs <- function(t, x, params) {
    A <-
      xtabs(interaction.coef ~ i + j, data = alpha.df[alpha.df$time == ceiling(t), ])
    # A <- alpha * (unique(env)[ceiling(t)]*10)
    # print(paste(ceiling(t), A[2,1]))
    x[x < Nmin] = 0
    dx = x * (r + A %*% x) # add a constant here
    
    # dx[x<Nmin]=0 # extinct species don't grow anymore - remove to allow for reinvasions
    return(list(dx))
  }
  
  result = ode(y = No, times = time, func = eqs)
  comm_dynamic <-
    gather(as.data.frame(result),
           key = "i",
           value = "n",
           -time)
  comm_dynamic$i <- as.factor(comm_dynamic$i)
  
  alpha.df$i <- as.factor(alpha.df$i)
  cwm_dynamic <-
    left_join(comm_dynamic, unique(alpha.df[, c("i", "trait.i")])) %>%
    group_by(time) %>%
    summarize("cwm" = mean( trait.i   * (n/sum(n))),
              "cwv" = sum ( trait.i^2 * (n/sum(n))) - mean( trait.i   * n/sum(n)))
  cwm_dynamic$env <- unique(alpha.df$env)
  
  
  
  if (plot.species.dynamics == T) {
    require(gridExtra)
    env.plot <-  ggplot(cwm_dynamic) +
      geom_line(aes(x = time, y = env)) +
      scale_color_viridis_d() +
      theme_classic()
    
    comm.plot <- ggplot(comm_dynamic, aes(y = n, x = time)) +
      geom_line(aes(col = i)) +
      scale_color_viridis_d() +
      theme(legend.position = "bottom")+
      theme_classic()

  print(grid.arrange(env.plot, comm.plot))
  }
  
  
  if (plot.response.diagram == T) {
    require(gridExtra)
    cwm.plot <-  ggplot(cwm_dynamic) +
      geom_line(aes(x = time, y = cwm)) +
      scale_color_viridis_d() +
      theme_classic()
    
    response.diagram <-  ggplot(cwm_dynamic) +
      geom_path(aes(x = env, y = cwm, col = time)) +
      scale_color_viridis() +
      # theme()+
      theme_classic(legend.position = "bottom")
    
    print(grid.arrange(env.plot, cwm.plot, response.diagram))
  }
return(left_join(comm_dynamic, cwm_dynamic))
}
