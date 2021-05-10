niche_difference      <- function(x) {
  ifelse(x$trait.i != x$trait.j, 
        (-1 / (1 + x$env)) * abs (1 / (x$trait.i - x$trait.j)),
         -1)
}
