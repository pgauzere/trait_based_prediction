extreme_facilitation  <- function(x) {
  ifelse(x$trait.i != x$trait.j, 
         1/(1 + x$env) * abs(x$trait.i - x$trait.j) / 
           max(1/(1 + x$env) * abs(x$trait.i - x$trait.j)),#this is the function that matters
         -1)
}