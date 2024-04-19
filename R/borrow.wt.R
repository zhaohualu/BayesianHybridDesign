#' Borrowing weights from historical study
#'
#' This function calculates the borrowing weights for the hybrid control from historical study.
#' The following methods are included: (1) Empirical Bayes: dynamic weight is determined by maximizing marginal distribution;
#' (2) Bayesian P value: similarity measured by Bayesian P value comparing two distributions;
#' (3) Generalized Bhattacharyya Coefficient (BC): Bhattacharyya coefficient is a special case for theta = 0.5.
#' (4) Jensen-Shannon divergence (JSD): a measurement of similarity for two distributions.
#'
#' @param Yc A scalar number of responders for control arm in current study
#' @param nc A scalar number of patients in control arm in current study
#' @param Ych A scalar number of responders for control treatment in historical study
#' @param nch A scalar total number of patients in historical control
#' @param nche A scalar for the equivalent number of patients borrowed from historical study
#' @param a0c The 1st scale hyperprior hyperparameter for control response rate beta(a0c, b0c). Default 0.001.
#' @param b0c The 2nd hyperprior hyperparameter for control response rate beta(a0c, b0c). Default 0.001.
#' @param delta_threshold A scale threshold parameter. Only if abs(pc (current study) - pch) <= delta_threshold, we borrow from historical control. Default 0.1.
#' @param method A string characters. Method for dynamic borrowing, "Empirical Bayes", "Bayesian p", "Generalized BC", "JSD"
#' @param theta A scalar parameter with a range of (0, 1), and applicable to "Generalized BC". Default 0.5.
#' @param eta A scalar parameter with a range of (0, infty), and applicable to methods "Bayesian p", "Generalized BC", "JSD". Default 1.
#'
#' @return An object with values
#'  \itemize{
#'  \item a Global borrowing weight
#'  \item wd Dynamic borrowing weight according to similarly of response rate
#'  \item w Overall borrowing weight
#'  }
#' @examples
#'
#' borrow.wt(Yc=40*0.3, nc=40, Ych=180*0.3, nch=180, nche=40, a0c=0.001, b0c=0.001, delta_threshold=0.1, eta=1)
#'
#' @export
#'
borrow.wt = function (Yc, nc,
                      Ych, nch, nche,
                      a0c=0.001, b0c=0.001,
                      delta_threshold=0.1,
                      method="Empirical Bayes", theta=0.5, eta=1){

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
  }else if (method =="Empirical Bayes a") {
    #Empirical weight
    wd = optimize(logwfunction, c(0,1),
                  a0c = a0c,
                  b0c = b0c,
                  Yc = Yc,
                  nc = nc,
                  Ych = a*Ych,
                  nch = a*nch,
                  maximum = TRUE)$maximum
  } else if (method == "Bayesian p") {
    #Heterogeneity

    #P(p_c > p_ch)
    P_c_p_ch = function(y,ac,bc,ach,bch){
      p = pbeta(y,ac,bc,lower.tail=F)*dbeta(y,ach,bch)
      Idx = is.infinite(p)
      p[Idx] = sign(p[Idx]) * 1e10
      return(p)
    }
    #P(p_ch > p_c)
    P_ch_p_c = function(y,ac,bc,ach,bch){
      p=pbeta(y,ach,bch,lower.tail=F)*dbeta(y,ac,bc)
      Idx = is.infinite(p)
      p[Idx] = sign(p[Idx]) * 1e10
      return(p)
    }
    xi1 = integrate(P_c_p_ch, lower=0.0001, upper=0.9999,
                    ac = ac, bc = bc, ach = ach, bch=bch)$value
    xi2 = integrate(P_ch_p_c, lower=0.0001, upper=0.9999,
                    ac = ac, bc = bc, ach = ach, bch=bch)$value
    wd = (2*min(xi1, xi2))^eta
  } else if (method == "Generalized BC") {
    # \int_0^1 \sqrt{f_1(x)f_2(x)}dx
    f.den = function(y){
      L1 = exp(theta*log(dbeta(y,ach,bch))+(1-theta)*log(dbeta(y,ac,bc)))
      L2 = exp(theta*log(dbeta(y,ac,bc))+(1-theta)*log(dbeta(y,ach,bch)))
      L = (L1+L2)/2
      Idx = is.infinite(L)
      L[Idx] = sign(L[Idx]) * 1e10
      return(L)
    }
    wd = integrate(f.den, lower=0.0001, upper=0.9999)$value
    wd = wd^eta
  }else if (method == "JSD") {

    f.den = function(y){
      logfc = dbeta(y,ac,bc,log = T)
      logfch = dbeta(y,ach,bch,log=T)
      logfbar = log((exp(logfc) + exp(logfch))/2)


      ans = (logfc-logfbar)*exp(logfc) + (logfch-logfbar)*exp(logfch)
      # print("logfc")
      # print(logfc)
      # print("logfbar")
      # print(logfbar)
      # print("logfch")
      # print(logfch)
      # print(ans)

      Idx = is.infinite(ans)
      ans[Idx] = sign(ans[Idx]) * 1e10
      Idx = is.na(ans)
      ans[Idx] = mean(ans[!Idx])
      #
      return(ans)
    }

    wd = (1-0.5*integrate(f.den, lower=0.0001, upper=0.9999)$value)^eta

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
    print(method)
    stop("method not supported!")
  }

  #Overall weight
  w = wd * a * (abs(Yc/nc - Ych/nch) < delta_threshold)

  return(list(a=a, wd=wd, w=w))
}


