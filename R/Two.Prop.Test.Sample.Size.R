#' Sample Size Calibration for Two-Sample Proportions
#'
#' This function calculates the sample size for comparing two-sample proportions
#' https://www2.ccrb.cuhk.edu.hk/stat/proportion/Casagrande.htm
#'
#'p1,p2,alpha,beta,r
#' @param p1 Arm 1 proportion (control)
#' @param p2 Arm 2 proportion (experimental)
#' @param alpha Type I error rate (one-sided)
#' @param beta type II error rate (1-power)
#' @param r Ratio of Arm 2 (experimental) vs Arm 1 (control)
#'
#' @return sample size for arm 1 (control) and total
#' @examples
#'
#' Two.Prop.Test.Sample.Size (p1=0.27, p2=0.47, alpha = 0.10, beta=0.2, r=2)
#'
#' @export
#'

Two.Prop.Test.Sample.Size = function(p1,p2,alpha,beta,r){

  delta = abs(p1-p2)
  Pbar = (p1+ r*p2)/(r+1)
  Qbar = 1-Pbar

  mp = (qnorm(alpha,lower.tail = F)*sqrt((r+1)*Pbar*Qbar)+
          qnorm(beta,lower.tail = F)*sqrt(r*p1*(1-p1)+p2*(1-p2)))^2/(r*delta^2)
  m = mp * (1+sqrt(1+2*(r+1)/(r*mp*delta)))^2/4

  return(c(ceiling(m),(r+1)*ceiling(m)))
}


