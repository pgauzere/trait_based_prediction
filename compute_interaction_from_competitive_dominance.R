competitive_dominance <- function(x) {
  ifelse(x$trait.i > x$trait.j,
         (1/(1 + x$env)) * (-(x$trait.i - x$trait.j)) /
           abs(min( (1/(1 + x$env)) * (-(x$trait.i - x$trait.j)))),
         ifelse(x$trait.i < x$trait.j, 0, -1 ))
}
