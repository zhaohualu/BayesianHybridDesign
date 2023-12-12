#' Borrowing weights from historical study
#' 
#' This function calculates the borrowing weights for the hybrid control from historical study
#' 
#' @param Yc number of responders for control arm in current study
#' @param nc number of patients in control arm in current study
#' @param Ych number of responders for control treatment in historical study
#' @param nch Total number of patients in historical control
#' @param nche Equivalent number of patients borrowed from historical study
#' @param a0c hyperprior for control response rate beta(a0c, b0c)
#' @param b0c hyperprior for control response rate beta(a0c, b0c)
#' @param delta_threshold Borrow when abs(pc (current study) - pch) <= delta_threshold
#'  
#' @return An object with values
#'  \itemize{
#'  \item a Global borrowing weight
#'  \item we Dynamic borrowingweight according to similarly of response rate
#'  \item w Overall borrowing weight
#'  }
#' @examples
#' 
#' borrow.wt(Yc=12, nc=40, Ych=70, nch=200, nche=40, a0c=0.001, b0c=0.001, delta_threshold=0.1)
#' 
#' @export
#' 
borrow.wt = function (Yc=12, nc=40,
                      Ych=70, nch=200, nche=40, 
                      a0c=0.001, b0c=0.001,
                      delta_threshold=0.1){

  #Target function to optimize weight
  logwfunction = function(w,a0,b0,Y,n,Y0,n0){
    #a0, b0: beta hyper prior parameters beta(a0, b0)
    #Y, n: current study number of responders and size in control
    #Y0, n0: historical control for number of responders and size
    
    lbeta(a0 + Y + w*Y0, b0 + (n-Y) + w*(n0-Y0)) - 
      lbeta(a0  + w*Y0, b0 + w*(n0-Y0))
  }
  
  #Empirical weight
  we = optimize(logwfunction, c(0,1),
                a0 = a0c,
                b0 = b0c,
                Y = Yc,
                n = nc, 
                Y0 = Ych,
                n0 = nch,
                maximum = TRUE)$maximum
  
  #Global weight parameter for borrowing    
  a = nche/nch
  
  #Overall weight
  w = we * a * (abs(Yc/nc - Ych/nch) < delta_threshold)
  
  return(list(a=a, we=we, w=w))
}


