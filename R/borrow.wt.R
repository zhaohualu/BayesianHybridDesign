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
#' @param method Method for dynamic borrowing, "Empirical Bayes" or "Heterogeneity".
#'   
#' @return An object with values
#'  \itemize{
#'  \item a Global borrowing weight
#'  \item we Dynamic borrowingweight according to similarly of response rate
#'  \item w Overall borrowing weight
#'  }
#' @examples
#' 
#' borrow.wt(Yc=40*0.312, nc=40, Ych=234*0.312, nch=234, nche=40, a0c=0.001, b0c=0.001, delta_threshold=0.1)
#' 
#' @export
#' 
borrow.wt = function (Yc=40*0.312, nc=40,
                      Ych=234*0.312, nch=234, nche=40, 
                      a0c=0.001, b0c=0.001,
                      delta_threshold=0.1, method="Empirical Bayes"){

  #Global weight parameter for borrowing    
  a = nche/nch
  
  #Target function to optimize weight
  logwfunction = function(w,a0,b0,Y,n,Y0,n0){
    #a0, b0: beta hyper prior parameters beta(a0, b0)
    #Y, n: current study number of responders and size in control
    #Y0, n0: historical control for number of responders and size
    
    lbeta(a0 + Y + w*Y0, b0 + (n-Y) + w*(n0-Y0)) - 
      lbeta(a0  + w*Y0, b0 + w*(n0-Y0))
  }
  
  if (method =="Empirical Bayes") {
    #Empirical weight
    we = optimize(logwfunction, c(0,1),
                  a0 = a0c,
                  b0 = b0c,
                  Y = Yc,
                  n = nc, 
                  Y0 = Ych,
                  n0 = nch,
                  maximum = TRUE)$maximum
  } else if (method == "Heterogeneity") {
    #Heterogeneity
    
    #P(p_c > p_ch)
    P_c_p_ch = function(y,ac,bc,ach,bch){
      pbeta(y,ac,bc,lower.tail=F)*dbeta(y,ach,bch)
    }
    #P(p_ch > p_c)
    P_ch_p_c = function(y,ac,bc,ach,bch){
      pbeta(y,ach,bch,lower.tail=F)*dbeta(y,ac,bc)
    }
    xi1 = integrate(P_c_p_ch, lower=0.0001, upper=0.9999,
                    ac = a0c + Yc, bc = b0c + (nc - Yc),
                    ach = a0c + a*Ych, bch=b0c + a*(nch - Ych))$value
    xi2 = integrate(P_ch_p_c, lower=0.0001, upper=0.9999,
                    ac = a0c + Yc, bc = b0c + (nc - Yc),
                    ach = a0c + a*Ych, bch=b0c + a*(nch - Ych))$value
    we = 2*min(xi1, xi2)
  }
  
  #Overall weight
  w = we * a * (abs(Yc/nc - Ych/nch) < delta_threshold)
  
  return(list(a=a, we=we, w=w))
}


