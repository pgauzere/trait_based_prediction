niche_difference <- function(x) {
trait_diff <- abs(x$trait.i - x$trait.j)

threshold <- quantile(trait_diff, probs = 0.1)

trait_diff[trait_diff < threshold] <- threshold

ifelse(x$trait.i != x$trait.j,
       -1/((1 + x$env) * trait_diff) / #this is the function that matters
         abs(min(-1/((1 + x$env) * trait_diff))),
       -1)
}
