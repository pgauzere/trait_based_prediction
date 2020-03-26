niche_difference      <- function(x) {
  ifelse(x$trait.i != x$trait.j, 
         -1/(1 + x$env) * abs(x$trait.i - x$trait.j) / #this is the function that matters
           abs(min(-1/(1 + x$env) * abs(x$trait.i - x$trait.j))),
         -1)
}