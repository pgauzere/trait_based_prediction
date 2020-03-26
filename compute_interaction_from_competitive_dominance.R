competitive_dominance <- function(x) {
  ifelse(x$trait.i > x$trait.j, 
         (1 + x$env) * -abs(x$trait.i - x$trait.j)/ #this is the function that matters
           max((1 + x$env) * (x$trait.i - x$trait.j)), 
         ifelse(x$trait.i < x$trait.j, 0, -1 ))
}
