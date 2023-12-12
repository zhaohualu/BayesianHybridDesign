#' Power Calculation by Fisher's Exact Test
#' 
#' This function calculates the power by simulations using Fisher's exact test.
#' 
#' @param pt Probablity of success in experimental arm
#' @param pc Probablity of success in control arm
#' @param nt Number of subject in experimental arm
#' @param nc Number of subject in control arm
#' @param alpha Type I error (one-sided)
#' @param nsim Number of replications in simulation
#' @param seed seed for simulation
#'  
#' @return power
#' 
#' @examples
#' 
#' fisher.power(pt = 0.512, nt = 40, pc = 0.312, nc = 40,nsim = 1000)
#' 
#' @export
#' 
fisher.power <- function(pt = 0.512, nt = 40, pc = 0.312, nc = 40, alpha=0.05,
                         nsim = 10000, seed=2023){
  set.seed(seed)
  sig=0
  for (i in 1:nsim){
    Yc = rbinom(1, size=nc, prob=pc)
    Yt = rbinom(1, size=nt, prob=pt)
    p = fisher(Yc=Yc, nc=nc, Yt=Yt, nt=nt, 
               alternative = "greater")$p.value
    if (p <= alpha) {sig = sig + 1}    
  }
  power = sig/nsim
  return(power)
}
