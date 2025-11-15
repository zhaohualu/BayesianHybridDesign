#' Statistical Analysis for a Bayesian Hybrid Design
#'
#' This function performs the statistical analysis for a Bayesian Hybrid Design using dynamic power prior approach.
#'
#' @param Yt Number of responses in experimental arm in current study
#' @param nt number of patients in experimental arm in current study
#' @param Yc Number of responses in control arm in current study
#' @param nc number of patients in control arm in current study
#' @param nche Equivalent number of patients borrowed from historical study
#' @param nch Total number of patients in historical control
#' @param Ych Number of responses in control treatment in historical study
#' @param sig significance boundary defined as P(pt_hat > pc_hat|hybrid) > sig
#' @param credlev: Credible interval level
#' @param a0c hyperprior for control response rate beta(a0c, b0c)
#' @param b0c hyperprior for control response rate beta(a0c, b0c)
#' @param a0t hyperprior for experimental response rate beta(a0t, b0t)
#' @param b0t hyperprior for experimental response rate beta(a0t, b0t)
#' @param delta_threshold Borrow when abs(pc_hat (current study) - pch) <= delta_threshold
#'
#' @return An object with values
#'  \itemize{
#'  \item prob.pt.gt.pc Probability of experimental arm having a better posterior response rate than control
#'  \item median_hca Posterior median response rate for hybrid control
#'  \item CI_hca Credible interval for median response rate for hybrid control
#'  \item median_c Posterior median response rate for current study control
#'  \item CI_c Credible interval for median response rate for current study control
#'  \item median_t Posterior median response rate for current study experimental arm
#'  \item CI_t Credible interval for median response rate for current experimental arm
#'  \item delta.m Posterior median response rate difference for experimental arm vs control arm based on hybrid design
#'  \item delta.CI Credible interval for posterior median response rate difference for experimental arm vs control arm  based on hybrid design
#'  \item delta.m_trial Posterior median response rate difference for experimental arm vs control arm based on current study only
#'  \item delta.CI_trial  Credible interval for posterior median response rate difference for experimental arm vs control arm based on current study only
#'  \item conclusion Statistical inference conclusion for comparing the response rate of experimental arm vs control arm
#'  }
#' @examples
#'
#' Bayesian.Hybrid.Analysis(Yt=18, nt=40,Yc=13,nc=40,Ych=73,nche=40,nch=234)
#'
#' @export
#'
Bayesian.Hybrid.Analysis = function(Yt=20, nt=40,Yc=12,nc=40,Ych=73,
                                    nche=40,nch=234, sig=0.90, credlev = 0.8,
                                    a0c=0.001,b0c=0.001,a0t=0.001,b0t=0.001,
                                    delta_threshold=0.1){

  #Dynamic borrowing parameter
  wt = borrow.wt(Yc=Yc, nc=nc, Ych=Ych, nch=nch, nche=nche,
                 a0c=a0c, b0c=b0c, delta_threshold=delta_threshold)
  w = wt$w

  #Poster distributions
  apost_c_trial = a0c + Yc
  bpost_c_trial = b0c + (nc-Yc)

  apost_c_hca = apost_c_trial + w*Ych
  bpost_c_hca = bpost_c_trial + w*(nch-Ych)

  apost_t = a0t + Yt
  bpost_t = b0t + (nt-Yt)

  P_ORRt_upper_Times_p_ORRc = function(y,ac,bc,at,bt){
    pbeta(y,at,bt,lower.tail=F)*dbeta(y,ac,bc)
  }
  phat_pt_larger_pc=integrate(P_ORRt_upper_Times_p_ORRc,
                              lower=0.0001, upper=0.9999,
                              ac = apost_c_hca, bc = bpost_c_hca,
                              at=apost_t,bt=bpost_t)$value
  if (phat_pt_larger_pc >= sig) {
    conclusion = "Statistically significant"} else {
      conclusion = "Not Statistically significant"
    }

  median_hca = qbeta(0.5, apost_c_hca, bpost_c_hca)
  CI_hca = c(qbeta(0.1, apost_c_hca, bpost_c_hca),
             qbeta(0.9, apost_c_hca, bpost_c_hca))
  median_c = qbeta(0.5, apost_c_trial, bpost_c_trial)
  CI_c = c(qbeta(0.1, apost_c_trial, bpost_c_trial),
           qbeta(0.9, apost_c_trial, bpost_c_trial))
  median_t = qbeta(0.5, apost_t, bpost_t)
  CI_t = c(qbeta(0.1, apost_t, bpost_t),
           qbeta(0.9, apost_t, bpost_t))

  #Calculate the posterior median of (pt_hat - pc_hat)
  #P(pt_hat - pc_hat < m) = 0.5
  mfun = function(m){
    yfun = function(y,ac,bc,at,bt){
      pbeta(m+y,at,bt) * dbeta(y,ac,bc)
    }
    ans = integrate(yfun,
                    lower=0, upper=1,
                    ac = apost_c_hca, bc = bpost_c_hca,
                    at=apost_t, bt=bpost_t)$value - 0.5
    return(ans)
  }

  delta.m = uniroot(f=mfun, interval = c(0, 1))$root

  #Calculate the posterior median of (pt_hat - pc_hat) for trial data only
  #P(pt_hat - pc_hat < m) = 0.5
  mfun_trial = function(m){
    yfun = function(y,ac,bc,at,bt){
      pbeta(m+y,at,bt) * dbeta(y,ac,bc)
    }
    ans = integrate(yfun,
                    lower=0, upper=1,
                    ac = apost_c_trial, bc = bpost_c_trial,
                    at=apost_t, bt=bpost_t)$value - 0.5
    return(ans)
  }

  delta.m_trial = uniroot(f=mfun_trial, interval = c(0, 1))$root

  #Calculate the credible interval of delta = (pt_hat - pc_hat)
  #P(delta.m - d < pt_hat-pc_hat < delta.m + d) = credlev
  #Let x = pt_hat, y = pc_hat, representing two random variables
  #P(delta.m - d < x - y < delta.m + d) = P(delta.m-d+y < x <delta.m + d+y)

  dfun = function(d){
    yfun = function(y,ac,bc,at,bt){
      (pbeta(delta.m+y+d,at,bt) - pbeta(delta.m+y-d,at,bt) )*dbeta(y,ac,bc)
    }
    ans = integrate(yfun,
                    lower=0, upper=1,
                    ac = apost_c_hca, bc = bpost_c_hca,
                    at=apost_t, bt=bpost_t)$value - credlev
    return(ans)
  }
  delta.d = uniroot(f=dfun, interval = c(0, 1))$root

  delta.CI = c(delta.m - delta.d, delta.m + delta.d)

  #Calculate the credible interval of (pt_hat - pc_hat) trial data only

  dfun_trial = function(d){
    yfun = function(y,ac,bc,at,bt){
      (pbeta(delta.m+y+d,at,bt) - pbeta(delta.m+y-d,at,bt) )*dbeta(y,ac,bc)
    }
    ans = integrate(yfun,
                    lower=0, upper=1,
                    ac = apost_c_trial, bc = bpost_c_trial,
                    at=apost_t, bt=bpost_t)$value - credlev
    return(ans)
  }
  delta.d_trial = uniroot(f=dfun_trial, interval = c(0, 1))$root

  delta.CI_trial = c(delta.m_trial - delta.d_trial, delta.m_trial + delta.d_trial)

  return(list(prob.pt.gt.pc=phat_pt_larger_pc, conclusion = conclusion,
              median_hca=median_hca, CI_hca=CI_hca,
              median_c=median_c, CI_c=CI_c,
              median_t=median_t, CI_t=CI_t,
              delta.m=delta.m, delta.CI=delta.CI,
              delta.m_trial=delta.m_trial, delta.CI_trial=delta.CI_trial))
}
