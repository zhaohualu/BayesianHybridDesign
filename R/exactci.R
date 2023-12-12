#' Clopper-Pearson Exact Confidence Interval
#' 
#' This function calculates the Clopper-Pearson exact confidence interval 
#' 
#' @param r Number of reponses
#' @param n Total number of subjects
#' @param conflev Confidence level. Always two-sided.
#'  
#' @return Clopper-Pearson CI
#' @examples
#' 
#' exactci (r=4, n=20, conflev=0.95)
#' 
#' @export
#' 
exactci <- function(r, n, conflev){
  alpha = (1 - conflev)
  if (r == 0) {
    ll = 0
    ul = 1 - (alpha/2)^(1/n)
  }
  else if (r == n) {
    ll = (alpha/2)^(1/n)
    ul = 1
  }
  else {
    ll = 1/(1 + (n - r + 1) / (r * qf(alpha/2, 2 * r, 2 * (n-r+1))))
    ul = 1/(1 + (n - r) / ((r + 1) * qf(1-alpha/2, 2 * (r+1), 2 *(n-r))))
  }
  c(ll,ul)
}
