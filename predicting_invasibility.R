
predict_invasibility_surface <- function(
  nsp = 10,
  env = seq(from = 0, to = 10,  by =1 ),  
  trait.distribution = "uniform",  
  mechanism = "competitive dominance" ,
  Nmin = 0.0001 ,
  initial.abundance = 0.0005 ,
  growth.rate = 0.5,
  n.time.step = 250 ,
  extinction = F ,
  plot.invasion = F,
  traits_to_test = seq(0,1,0.1)){



nsp = 10
env = seq(from = 0, to = 10,  by = 1)
trait.distribution = "uniform"
mechanism = "niche difference"
Nmin = 0.0001
initial.abundance = 0.0005
growth.rate = 0.5
n.time.step = 250
extinction = T
plot.invasion = T
traits_to_test = seq(0,1, 0.1)


sub_alpha <- predict_demographic_model_parameters(rep = 1, 
                                                 nsp = nsp, 
                                                 env = env, 
                                                 trait.distribution = trait.distribution, 
                                                 mechanism = mechanism)


nsp <- length(unique(sub_alpha$i)) + 1
env <- unique(sub_alpha$env)

# create empty matrices
alpha <- array(rep(NA, nsp^2), dim = c(nsp, nsp, length(env), length(mechanism), length(trait.distribution)),  
               dimnames = list(1:nsp, 1:nsp, env, mechanism, trait.distribution))

# reshape matrix to dataframe (for computation purpose)
alpha.df2 <- reshape2::melt(alpha)
colnames(alpha.df2) <- c("i", "j", "env", "mechanism", "trait.distribution", "interaction.coef")
alpha.df2$trait.i <- NA
alpha.df2$trait.j <- NA
alpha.df2[alpha.df2$i !=11 & alpha.df2$j !=11, c("interaction.coef", "trait.i", "trait.j")] <- sub_alpha[, c("interaction.coef", "trait.i", "trait.j")]
alpha.df2[alpha.df2$i !=11 & alpha.df2$j ==11, "trait.i"] <- rep(unique(sub_alpha$trait.i))
alpha.df2[alpha.df2$i ==11 & alpha.df2$j !=11, "trait.j"] <- rep(unique(sub_alpha$trait.j))


res_invasion <-
  matrix(
    data = NA,
    nrow = length(env),
    ncol = length(traits_to_test),
    dimnames = list(as.character(env),
                    as.character(traits_to_test))
  )


for(i in 1:length(env)){

  # alpha.stability <- plyr::ddply(alpha.df, .(mechanism, trait.distribution, replicate), function(x){
  x <- sub_alpha[sub_alpha$env == env[i],]
  
  ### Fist simulation
  # initial abundance
  if (class(initial.abundance) == "numeric") {
    if (length(initial.abundance) == 1) {
      No = rep(initial.abundance, nsp-1)
    } else{
      No = initial.abundance
    }
  } else {
    if (initial.abundance == "random") {
      No = runif(nsp-1)
    }
  }
  
  # growth rate
  if (class(growth.rate) == "numeric") {
    if (length(growth.rate) == 1) {
      r = rep(growth.rate, nsp-1)
    } else{
      r = growth.rate
    }
  } else {
    if (growth.rate == "random") {
      # intrinsic growth rates sampled from normal distribution
      r = rnorm(nsp-1, 1, 0.1)
    }
    if (growth.rate == "trait") {
      r = x$trait.i[1:nsp-1]
    }
  }
  
  r_base<-r  
  
  #define equation system
  eqs <- function(t,x,params){
    x[ x < Nmin ] = 0 
    dx = x*(r + A%*%x ) # add a constant here 
    dx[x < Nmin] = 0 # extinct species don't grow anymore - remove to allow for reinvasions
    return(list(dx))
  }  
  
  ### Assess stability of A matrix  
  #convert table into matrix
  A <- xtabs(interaction.coef ~ i + j, data=x)
  
  ### Numerical simulations
  #run ODE dynamics
  result = ode(y = No, times=1:250, func=eqs)
  comm_final <- result[100,-1]
  
  for(j in 1:length(traits_to_test)){

    alpha.df2[alpha.df2$i == nsp, "trait.i"] <- traits_to_test[j]
    alpha.df2[alpha.df2$j == nsp, "trait.j"] <- traits_to_test[j]
    
    

    #compute interaction coeffcients depending on mechanisms and environment value
    alpha.df2$interaction.coef[(alpha.df2$i == nsp  |
                                  alpha.df2$j == nsp ) &
                                 (alpha.df2$mechanism == 'niche difference')] <-
      niche_difference(x = alpha.df2[(alpha.df2$i == nsp  |
                                        alpha.df2$j == nsp) &
                                       (alpha.df2$mechanism == 'niche difference'),])
    
    alpha.df2$interaction.coef[(alpha.df2$i == nsp  |
                                  alpha.df2$j == nsp) &
                                 (alpha.df2$mechanism == 'competitive dominance')] <-
      competitive_dominance(x = alpha.df2[(alpha.df2$i == nsp  |
                                             alpha.df2$j == nsp) &
                                            (alpha.df2$mechanism == 'competitive dominance'), ])
    
    alpha.df2$interaction.coef[(alpha.df2$i == nsp  |
                                  alpha.df2$j == nsp) &
                                 (alpha.df2$mechanism == 'extreme facilitation')] <-
      extreme_facilitation(x = alpha.df2[(alpha.df2$i == nsp  |
                                            alpha.df2$j == nsp) &
                                           (alpha.df2$mechanism == 'extreme facilitation'), ])
    
    
    
    
    
    x2 <- alpha.df2[alpha.df2$env==env[i] & alpha.df2$mechanism==mechanism & alpha.df2$trait.distribution==trait.distribution, ]
    
    No = c(comm_final, Nmin*2) # initial conditions
    #convert table into matrix
    A <- xtabs(interaction.coef ~ i + j, data=x2)
    # r <- c(r_base, runif(1,0,1))
    r <- c(r_base, traits_to_test[j])
    result = ode(y = No, times=1:n.time.step, func=eqs)
    
    # plot community dynamic
    comm_dynamic <- gather(as.data.frame(result), key="species_i", value="n", -time)

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
          labs(title = paste("trait value =",traits_to_test[j], "| env value =",env[i]))+
          ylim(0,1)
      )
    }
    
    res_invasion[i,j] <- comm_dynamic[nrow(comm_dynamic), "n"]

  #monitor progress
  print(paste("invasion : ","env[", i,"]/",length(unique(alpha.df$env)), "; trait[", j, "]/", length(traits_to_test),  "-- done", sep =''))
  
  }
}
  return(res_invasion)
}