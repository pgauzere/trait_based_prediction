require(reshape2)
create_empty_alpha.df <- function(nsp, env, mechanism, trait.distribution, rep){
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
  
  return(alpha.df)}