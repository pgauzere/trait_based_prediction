################################################################################
######################### Simulate Environmental Change ########################
########################################### Gauzere P. #########################

#This function simulate temporal environmental dynamics of different kinds
#It follow the basic principles proposed by Ryo et al., 2019 tree
#ref:
#Ryo, Masahiro, Carlos A. Aguilar-Trigueros, Liliana Pinek, Ludo AH Muller, and Matthias C. Rillig. "Basic principles of temporal dynamics." Trends in ecology & evolution (2019).

# function argument are : 
  # type [character] : "pulse", "step", "ramp", "linear, "cyclic", "non-stationary"
  # trend_value [numeric] : the slope of linear change, only works with type = "linear"
  # stochasticity_value [numeric] : the amount noise in the timeseries, only works with type = "linear", "cyclic", and "non-stationary"
  # cycle_value [numeric] : the amplitude of cyclic fluctuations, only works with type = "cyclic"
  # cycle_period [numeric] : the frequency of cyclic fluctuations, only works with type = "cyclic"
  # cycle_period [numeric] : the frequency of cyclic fluctuations, only works with type = "cyclic"
  # dispersion_value [numeric] : the amplitude of dispersion at each timesteps, only works with type = "non-stationary"

#function returns a dataframe with columns :
 # time [numeric] : the time steps
 # env.dynamic [numeric] : the environmental value at each time step

simulate_environmental_change <- function(type = "pulse",
                                          n.time.step = 100,
                                                  trend_value = 2,
                                                  stochasticity_value = 0,
                                                  cycle_value = 5,
                                                  cycle_period = 0.1, 
                                                  dispersion_value = 0.300,
                                                  plot.env.change = T) {
  
  
   time <- 1:n.time.step
 
  ##### simulate different temporal forcings
  
  ###discrete events
  ## pulse
  if (type == "pulse") {
    env.dynamic <- (exp( - ((time/(max(time)*0.05))-5)^2))*10
  }
  
  ## step
  if (type == "step") {
    env.dynamic <- (-10 / (1 + exp((time/(max(time)*0.05)) -10))) + 10 
  }
  
  ## ramp
  if (type == "ramp") {
    env.dynamic <- time*2
    env.dynamic[0:round(max(time)/3)] <- min(env.dynamic[round(max(time)/3) : max(time)])
    env.dynamic <- env.dynamic - min(env.dynamic)
    env.dynamic <- (env.dynamic *10) / max(env.dynamic)
  }
  
  ###trajectory
  if (type == "linear") {
    ## Trend stationary 
    env.dynamic <- scale(time) * trend_value + rnorm(length(time),sd=stochasticity_value) 
    env.dynamic <- env.dynamic + abs(min(env.dynamic))
  }
  
  ## Cyclo stationary 
  if (type == "cyclic") {
    env.dynamic <-  sin(time/(max(time)*cycle_period)) * cycle_value  + rnorm(length(time),sd=stochasticity_value)
    env.dynamic <- env.dynamic + abs(min(env.dynamic))
    plot(time, env.dynamic, type = "l" )
  }
  
  ## Non stationary 
  if (type == "non-stationary") {
    env.dynamic <-  cumsum(sample(c(-dispersion_value, dispersion_value), length(time), TRUE)) + rnorm(length(time),sd=stochasticity_value)
    env.dynamic <- env.dynamic + abs(min(env.dynamic))
  }
   
  environmental.change <- data.frame(time, env.dynamic)
  
  if(plot.env.change == T){
    require(ggplot2)
    print(
      ggplot(environmental.change, aes(time, env.dynamic)) +
      geom_line() +
      theme_classic()
      )}
  return(environmental.change)
    }
  