BC.CI <-
function(theta,sims,conf.level=0.95) 
{
  low <- (1 - conf.level)/2
  high <- 1 - low
  z.inv <- length(theta[theta < mean(theta)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(theta) - theta)
  top <- sum(U^3)
  under <- (1/6) * (sum(U^2))^{3/2}
  a <- top/under
  lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower2 <- lower <- quantile(theta, lower.inv)
  upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper2 <- upper <- quantile(theta, upper.inv)
  return(c(lower, upper))
}
