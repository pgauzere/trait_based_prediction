################################################################################
############ Predict feasfeacoexistence outcomes and network properties ###############
########################################### Gauzere P. #########################

# This is the whole enchilada. The function predict several indices of coexistence 
# and network properties from a given matrix (in dataframe long format) computed by predict_model_parameter() 

# alpha.df created by create_alpha_matrices.R has one community interaction matrix per replicat, env, trait distribution, mechanism
# an efficient way to use compute_alpha_properties is tidyverse on alpha.df dataframe
# e.g : 
# alpha.properties <-
#   alpha.df %>%
#   group_by(env, mechanism) %>%
#   group_modify( ~ compute_alpha_properties(.x))


# function requires a bunch of libraries to be installed : tidyverse, deSolve, rootSolve, igraph, ggraph, NetIndices

# function argument are : 
#  x [tibble]     : a subset of alpha.df created by A_create_alpha_matrices.R 
#  Nmin [numeric] : the relative abundance extinction threshold under which a population is extinct [numeric], between 0 and 1 
#  initial.abundance [multi] : the intial relative abundance of species. Can be : 
#            "random" [character] : randomly sampled from uniform distribution); 
#             a numeric value [numeric] : all species have the same initial abundance); 
#             a numeric vector [numeric] of length nsp (initial abundance defined for each species) 
#  growth.rate [multi] intrinsic growth rate of species. Can be : 
#             "random" [character] randomly sampled from normal distribution) ;
#             "trait" [character] growth rate of species is determined by its trait value;
#             a numeric value  [numeric] all species have the same intrinsic growth rate; 
#             a numeric vector of length nsp  [numeric] : growth rate defined for each species
#  extinction [logical]   : if TRUE, extinct species (i.e N < Nmin) don't grow anymore, if FALSE reinvasions are possible
#  n.time.step [numeric]  : the number of timse step for numeric simulations
#  plot.dynamic [logical] : if TRUE, plot the community dynamic
#  plot.network [logical] : if TRUE, plot the community network
#  network.threshold [numeric] : threshold value for interaction strenght under which the interaction is considered non-existent

# function returns a dataframe with columns : 
# "env" : environemtal value
# "mechanism" : assembly mechanism
#  "analytic_stability" = analytic_stability, based on the rule of "all negative eigenvalues" 
#  analytic_stability2" = analytic_stability2, based on the rule of "negative rightmost eigenvalue
#  "numeric_stability" = numeric_stability, based on variation of abundance at the end of numeric simulation
#  "analytic_feasibility" = analytic_feasibility,based on species abundances at analytic equilibrium 
#  "numeric_feasibility" = numeric_feasibility,based on species abundances at the end of numeric simulation
# "E1" : first moment of the distribution for the off-diagonal elements of the interaction matrix (i.e approximate the mean)
# "E2" : second distribution for the off-diagonal elements of the interaction matrix (i.e variance)
# "Ec" : third distribution for the off-diagonal elements of the interaction matrix (i,e ??)
          #E1, E2, Ec are supposed to predict the stability, see Grilli et al., Nature Comm 2017
# "connectance" : connectance of the netwrok 
# "modularity" : modularity
# "transitivity" : transitivity
# "n_realized"  : number of species with n > Nmin at t = n.time.steps
# "mean_trait" : the mean value of traits sampled to create the pool  
# "var_trait" : the variance of trait values sampled to create the pool
# "cwm" : the community abundance weighted mean trait value
# "cwv": : the community abundance weighted variance of trait value


predict_coexistence_outcome <- function(x,
                                     Nmin = 0.001,
                                     initial.abundance = "random",
                                     growth.rate = "trait",
                                     n.time.step = 250,
                                     extinction = F,
                                     plot.dynamic = F,
                                     plot.network = F, 
                                     network.threshold = 0.05) {
  
  # alpha.df <- predict_demographic_model_parameters(
  #   rep = 1,
  #   nsp = 20,
  #   env = seq(from = 0, to = 10, by = 0.1))
  # 
  # x <- alpha.df
  # Nmin = 0.001
  # initial.abundance = 0.01
  # growth.rate = 0.5
  # n.time.step = 250
  # extinction = F
  # plot.dynamic = F
  # plot.network = F
  # network.threshold = 0.05

### load library
  require(tidyverse)
  require(deSolve)
  require(rootSolve)
  require(igraph)
  require(ggraph)
  require(NetIndices)
  require(viridis)

### define ODE equations system for multi-species dynamics
  ##define parameters for ODE system
    #number of species (from alpha.df)
  nsp = length(unique(x$i))
  
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
  eqs <- function(t, x, params) {
    x[x < Nmin] = 0
    dx = x * (r + A %*% x) # add a constant here
    if(extinction == T){dx[x<Nmin]=0} else{} # extinct species don't grow anymore - remove to allow for reinvasions
    return(list(dx))
  }  
  

#### 1 Analytical study of dynamic system ####
  
  #### 1.1 Assess  stability of A matrix ####
  
  #convert table into matrix
  A <- xtabs(interaction.coef ~ i + j, data=x)
  
  #compute eigenvalues of Jacobian matrix 
  eigenJacobianA <- eigen(jacobian.full(y = No, func=eqs))$values
  
  #assess stability based on "all negative eigenvalues" and store  
  analytic_stability <- if (all(Re(eigenJacobianA) < 0)) {
    "stable"
  } else {
    "unstable"
  }
  
  #assess stability based on "negative rightmost eigenvalue of diag(xstar)*A" and store  
  #compute eigenvalues of A 
  analytic_stability2 <- if (min(Re(eigenJacobianA)) < 0) {
    "stable"
  } else {
    "unstable"
  }
  
  #### 1.2 Analyticaly assess feasability of A matrix ####
  
  # analytically assess species densities at equilibirum point  
  n_equilibrium_analytic = -solve(A) %*% r
  
  #assess feasibility based on "all species have abundance > 0 at equilibrium point" and store  
  analytic_feasibility <- if (all(n_equilibrium_analytic > Nmin)) {
    # all(n_equilibrium == -solve(A) %*% r)) {
    "feasible"
  } else {
    "nonfeasible"
  }

#### 2 Numerical study of dynamic system ####
  
  #run ODE dynamics
  result = ode(y = No, times=1:n.time.step, func=eqs)
  
  #### 2.1 Numerically assess stability of A matrix ####
  # stability after n time steps     
  dndt_equilibrium_numerical  = result[n.time.step,-1] - result[n.time.step-1, -1]
  
  numeric_stability <- if (all(dndt_equilibrium_numerical < Nmin)) {
    "stable"
  } else {
    "unstable"
  }
  
  #### 2.2 Numerically assess feasability of A matrix ####
  
  #assess feasibility based on "all species have abundance > 0"
  # numerically assess species densities after 1000 time steps
  n_equilibrium_numerical  = result[n.time.step,-1]
  
  numeric_feasibility <- if (all(n_equilibrium_numerical > Nmin)) {
    # all(n_equilibrium == -solve(A) %*% r)) {
    "feasible"
  } else {
    "nonfeasible"
  }
  
  #### 2.3 Size of realized community
  n_realized <- length(which(n_equilibrium_analytic > 0))
  comm_dynamic <- gather(as.data.frame(result), key="species_i", value="n", -time)
  
  # plot community dynamic
  if(plot.dynamic == T){
  print(
    ggplot(comm_dynamic, aes(y=n, x=time, col=species_i))+
          geom_line()+
          scale_color_viridis_d(guide=F)+
          theme_classic()+
          geom_hline(yintercept = Nmin, linetype = 2, col="red")+
          labs(title = paste("replicat #", as.character(x$replicat[1]),"\n",
                             "analytic: ", analytic_feasibility,"- 1:",analytic_stability, " | 2:", analytic_stability2, "\n",
                             "numerical: ", numeric_feasibility,"-",numeric_stability, "\n"
                             ))
        )
  }
  
  
###### 3 Compute network properties ######  
  
  #### 3.1 Build, prune, simplify network ####
  
  ## create network
  # define nodes
  nodes <- as.data.frame(unique(x$i)) #Set nodes to be the unique values of "i", i.e. the the trait values
  colnames(nodes) <- "id" #Change column name to id
  
  # define edges and weights
  edges <- x[, c('i', 'j', "interaction.coef")] #from the x dataframe define a edges dataframe, ie what individuals are interacting, which are the trait values and interaction coefs
  edges$weight <-  edges$interaction.coef
  
  ## prune network : interaction coefficient under a given value are deleted from network 
  interSpeInteractions <- edges$weight[edges$i != edges$j]
  # threshold <- abs(mean(interSpeInteractions) + qnorm(network.threshold)*sd(interSpeInteractions)/sqrt(length(interSpeInteractions)))
  edges <- edges %>% filter(abs(weight) < network.threshold)
  # edges <- edges %>% filter(abs(weight) < threshold)
  
  ## assemble network
  network <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE) #Directed is whether the inx between two indiviuals is bidirectional or not
  
  ##simplify network
  # We have 'backwards' duplicates in this file, i.e. a interacts with b and b interacts with a. This is not necessary since our networks are undirected. To check this....is_simple checks whether a graph is simple. simplify removes the loop and/or multiple edges from a graph. If both remove.loops and remove.multiple are TRUE the function returns a simple graph.
  simple_network <- igraph::simplify(network, edge.attr.comb="mean")
  
  #### 3.2 compute and store network properties ####
  ## connectance 
  # needs "GenInd()" function from netIndices (requires an input of an adjacency matrix)
  network.adj <- get.adjacency(network,sparse=F)
  connectance <- GenInd(network.adj)$C  # This function measures connectance as L/(N*(N-1)) where L is links, and N is nodes (can also be calculated as L/(N^2)
  
  ##compute modularity 
  cluster_network <- cluster_walktrap(simple_network) #"This function tries to find densely connected subgraphs, also called communities in a graph via random walks. The idea is that short random walks tend to stay in the same community." 
  modularity <- modularity(cluster_network)
  
  ##compute transitivity
  # "transitivity" = clustering coefficient measures the probability that the adjacent nodes of a node are connected. 
  transitivity <- transitivity(network, type = "global")  # gives the clustering coefficient of the whole network
  transitivity(network, type="local") # gives the clustering coefficient of each node 
  
  ## plot netwok
  if (plot.network == T) {
    print(
      ggraph(simple_network, layout = "linear") +
        geom_edge_arc(aes(
          width = abs(weight),
          color = weight,
          alpha = ..index..
        )) +
        geom_node_point(size = 4) +
        scale_edge_width(range = c(0, 4), guide = F) +
        labs(edge_width = "strengh", edge_color = "coefficient") +
        theme_graph() +
        scale_edge_color_viridis(
          "interaction coefficient",
          option = "A",
          limits = c(-1, 1),
          direction = -1
        ) +
        scale_edge_alpha('edge direction', guide = 'edge_direction')
    )
  }
  
  diag(A)<-NA
  E1         = sum(A, na.rm = T) / nsp*(nsp-1)
  E2         = sum(A^2 - E1^2, na.rm = T) / nsp*(nsp-1)
  Ec         = sum(A*t(A) - (E1^2/E2^2), na.rm = T) / (nsp*(nsp-1)*E2^2)
  
    return(data.frame(
    "analytic_stability" = analytic_stability, 
    "analytic_stability2" = analytic_stability2,
    "numeric_stability" = numeric_stability,  
    "analytic_feasibility" = analytic_feasibility,
    "numeric_feasibility" = numeric_feasibility,
    "E1"         = E1[[1]],
    "E2"         = E2[[1]],
    "Ec"         = Ec[[1]],
    "connectance" = connectance, 
    "modularity"  = modularity, 
    "transitivity" = transitivity, 
    "n_realized"= n_realized, 
    "mean_trait" = mean(x$trait.i),
    "var_trait" = var(x$trait.i),
    "cwm" = mean( unique(x$trait.i)   * n_equilibrium_numerical / sum(n_equilibrium_numerical)),
    "cwv" = sum ( unique(x$trait.i)^2 * n_equilibrium_numerical / sum(n_equilibrium_numerical)) - 
               mean( unique(x$trait.i)   * n_equilibrium_numerical / sum(n_equilibrium_numerical))
    ))
  }







