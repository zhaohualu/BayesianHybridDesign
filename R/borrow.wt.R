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
#' @param method Method for dynamic borrowing, "Empirical Bayes", "Bayesian p", "Expected Bayes Factor", "JSD"
#' @param theta A parameter with a range of (0, 1), and applicable to method = "Expected Bayes Factor" or "JSD", where the weight is defined as E_c((fch/fc)^theta) for "Expected Bayes Factor" and (1-0.5*(KL(fc|fbar)+KL(fch|fbar)))^(1/theta), and fbar=(fc+fch)/2.
#'
#' @return An object with values
#'  \itemize{
#'  \item a Global borrowing weight
#'  \item wd Dynamic borrowing weight according to similarly of response rate
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
                      delta_threshold=0.1, method="Empirical Bayes", theta=0.5){

  #Global weight parameter for borrowing
  a = nche/nch

  #Target function to optimize weight for Empirical Bayes method
  logwfunction = function(w,a0c,b0c,Yc,nc,Ych,nch){
    #a0c, b0c: beta hyper prior parameters beta(a0c, b0c)
    #Yc, nc: current study number of responders and size in control
    #Ych, nch: historical control for number of responders and size
    
    lbeta(a0c + Yc + w*Ych, b0c + (nc-Yc) + w*(nch-Ych)) -
      lbeta(a0c  + w*Ych, b0c + w*(nch-Ych))
  }
  
  ac = a0c + Yc
  bc = b0c + (nc - Yc)
  ach = a0c + a*Ych
  bch = b0c + a*(nch - Ych)

  if (method =="Empirical Bayes") {
    #Empirical weight
    wd = optimize(logwfunction, c(0,1),
                  a0c = a0c,
                  b0c = b0c,
                  Yc = Yc,
                  nc = nc,
                  Ych = Ych,
                  nch = nch,
                  maximum = TRUE)$maximum
  } else if (method == "Bayesian p") {
    #Heterogeneity

    #P(p_c > p_ch)
    P_c_p_ch = function(y,ac,bc,ach,bch){
      pbeta(y,ac,bc,lower.tail=F)*dbeta(y,ach,bch)
    }
    #P(p_ch > p_c)
    P_ch_p_c = function(y,ac,bc,ach,bch){
      pbeta(y,ach,bch,lower.tail=F)*dbeta(y,ac,bc)
    }
    xi1 = integrate(P_c_p_ch, lower=0, upper=1,
                    ac = ac, bc = bc, ach = ach, bch=bch)$value
    xi2 = integrate(P_ch_p_c, lower=0, upper=1,
                    ac = ac, bc = bc, ach = ach, bch=bch)$value
    wd = 2*min(xi1, xi2)
  } else if (method == "Expected Bayes Factor") {
    # \int_0^1 \sqrt{f_1(x)f_2(x)}dx
    f.den = function(y){
      L1 = exp(theta*log(dbeta(y,ach,bch))+(1-theta)*log(dbeta(y,ac,bc)))
      L2 = exp(theta*log(dbeta(y,ac,bc))+(1-theta)*log(dbeta(y,ach,bch)))
      L = (L1+L2)/2
      return(L)
    }
    wd = integrate(f.den, lower=0.0001, upper=0.9999)$value
  }else if (method == "JSD") {
    
    f.den = function(y){
      fc = dbeta(y,ac,bc)
      fch = dbeta(y,ach,bch)
      fbar = (fc + fch)/2
      
      ans = log(fc/fbar)*fc + log(fch/fbar)*fch
      return(ans)
    }
    eps = 1/theta
    wd = (1-0.5*integrate(f.den, lower=0.0001, upper=0.999999)$value)^eps
    
  } else if (method == "Density Product 2") {
    # \int_0^1 \sqrt{f_1(x)f_2(x)}dx
    f.den = function(y){
      (dbeta(y,ac,bc)*dbeta(y,ach,bch))
    }
    fc.den = function(y){
      (dbeta(y,ac,bc)^2)
    }
    fch.den = function(y){
      (dbeta(y,ach,bch)^2)
    }
    wd = integrate(f.den, lower=0.0001, upper=0.9999)$value/
      sqrt(integrate(fc.den, lower=0.0001, upper=0.9999)$value * integrate(fch.den, lower=0.0001, upper=0.9999)$value)

  } else{
    stop("method not supported!")
  }

  #Overall weight
  w = wd * a * (abs(Yc/nc - Ych/nch) < delta_threshold)

  return(list(a=a, wd=wd, w=w))
}


