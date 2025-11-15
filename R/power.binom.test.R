#' Power calculation using exact binomial test for one sample problem
#'
#' This function performs Fishers exact test customized for comparing response rate for two arms.
#'
#' @param n Sample size
#' @param p Assumed proportion of experimental treatment
#' @param p0 Historical proportion as a benchmark
#' @param alpha Type I error
#' @param alternative Alternative hypothesis: "two.sided", "less", "greater".
#' @param nsim Number of simulated trials, default 100,000
#' @param seed Seed for simulation
#'
#' @return Power
#' @examples
#'
#' power.binom.test (n=110, p=0.4, p0=0.128)
#'
#' @export
#'
power.binom.test = function(n=20, p=0.5, p0=0.3, alpha=0.05,
                            alternative = c("two.sided", "less", "greater"),
                            nsim = 100000, seed= 2025){
  set.seed(seed)
  
  x = rbinom(n=nsim, size=n, prob=p)
  p.value = rep(NA, nsim)
  
  for (i in 1:nsim){
    t = binom.test(x=x[i], n=n, p = p0, alternative = alternative)
    p.value[i] =t$p.value 
  }
  power = sum(p.value <= alpha) / nsim
  return(power)
}

