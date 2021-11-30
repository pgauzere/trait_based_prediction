niche_difference_non_normalized <- function(x) {
  ifelse(x$trait.i != x$trait.j, 
         (1 / (1 + x$env)) * (-1 / (abs(x$trait.i - x$trait.j))),
         -1)
}


niche_difference <- function(x) {
  trait_diff <- abs(x$trait.i - x$trait.j)
  
  threshold <- quantile(trait_diff, probs = 0.1)
  
  trait_diff[trait_diff < threshold] <- threshold
  
  ifelse(x$trait.i != x$trait.j,
         (1 / (1 + x$env)) * (-1/trait_diff) / #this is the function that matters
           abs(min((1 / (1 + x$env)) * (-1/trait_diff))),
         -1)
}