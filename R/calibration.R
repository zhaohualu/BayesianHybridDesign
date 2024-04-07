#' Calibration for Bayesian Hybrid Design
#' 
#' This function calculates the threshold tau for calibration of Bayesian Hybrid Design.
#' P(pt>pc|hybrid data)>tau is used to determin statistical significance.
#' 
#' @param nt number of patients in experimental arm in current study    
#' @param nc number of patients in control arm in current study  
#' @param nche Equivalent number of patients borrowed from historical study
#' @param nch Total number of patients in historical control
#' @param alpha Type I error rate. The Bayesian decision threshold tau is calibrated towards this rate. 
#' @param pc Response rate for control arm in current study
#' @param pch Response rate for control treatment in historical study
#' @param a0c Hyperprior for control response rate beta(a0c, b0c)
#' @param b0c Hyperprior for control response rate beta(a0c, b0c)
#' @param a0t Hyperprior for experimental response rate beta(a0t, b0t)
#' @param b0t Hyperprior for experimental response rate beta(a0t, b0t)
#' @param delta_threshold Borrow when abs(pc_hat (current study) - pch) <= delta_threshold 
#' @param nsim Number of replications to calculate power
#' @param seed=2024 seed for simulations
#' @param method Method for dynamic borrowing, "Empirical Bayes", "Bayesian p", "Generalized BC", "JSD"
#' @param theta A parameter with a range of (0, 1), and applicable to "Generalized BC".
#' @param eta A parameter with a range of (0, infty), and applicable to methods "Bayesian p", "Generalized BC", "JSD". Default 1.
#' @param datamat A matrix with dimension nsim * 2 containing the pre-simulated data for the study treatment (1st column) and control  (1st column) groups, respectively. If not supplied, binomial random Monte Carlo samples will be generated in the function.
#' @param w0 prior power parameters w. If not specified (default), w_d is calculated by the specified method for dynamic borrowing.
#' 
#' @return The threshold for statistical significance that can control the type I error (1-sided)
#' @examples
#' 
#' calibration (nt=45, pc=0.312, nc=45, pch=0.312, nche=45, nch=234, alpha = 0.10, a0c=0.001, b0c=0.001, a0t=0.001, b0t=0.001, delta_threshold=0.1, method="Empirical Bayes", nsim = 100000, seed=2024)
#' calibration (nt=45, pc=0.312, nc=45, pch=0.312, nche=45, nch=234, alpha = 0.10, a0c=0.001, b0c=0.001, a0t=0.001, b0t=0.001, delta_threshold=0.1, method="Generalized BC", theta=0.5, eta=1, nsim = 100000, seed=2024)
#' calibration (nt=45, pc=0.312, nc=45, pch=0.312, nche=45, nch=234, alpha = 0.10, a0c=0.001, b0c=0.001, a0t=0.001, b0t=0.001, delta_threshold=0.1, method="Bayesian p", theta=NULL, eta=1, nsim = 100000, seed=2000)
#' calibration (nt=45, pc=0.312, nc=45, pch=0.312, nche=45, nch=234, alpha = 0.10, a0c=0.001, b0c=0.001, a0t=0.001, b0t=0.001, delta_threshold=0.1, method="JSD", theta=NULL, eta=2, nsim = 100000, seed=2000)
#' 
#' @export
#' 
calibration = function(nt=40, pc=0.312, nc=40,pch=0.312,
                       nche=40,nch=234, alpha = 0.10,
                       a0c=0.001,b0c=0.001,a0t=0.001,b0t=0.001,
                       delta_threshold=0.1, 
                       method="Empirical Bayes", theta=0.5,eta=1,
                       datamat = NULL, w0 = NULL,
                       nsim = 10000, seed=2000){
  
  #under H0
  pt = pc
    
  #Generate P(pt>pc|hybrid data) for each simulated trial
  set.seed(seed)
  pt_larger_pc = rep(NA, nsim)
  Ych=nch*pch
  
  for(i in 1:nsim){
    if(is.null(datamat)){
      Yt.s = rbinom(1,size=nt,prob=pt) #simulated Yt
      Yc.s = rbinom(1,size=nc,prob=pc) #simulated Yc
    } else{
      Yt.s =datamat[i,1]
      Yc.s =datamat[i,2]
    }

    #Dynamic borrowing parameter
    if(is.null(w0)){
      wt = borrow.wt(Yc=Yc.s, nc=nc, Ych=Ych, nch=nch, nche=nche,  
                     a0c=a0c, b0c=b0c, delta_threshold=delta_threshold, 
                     method=method, theta=theta, eta=eta)
      w = wt$w
    } else {
      w=w0
    }
    
    #Posterior distributions
    apost_c_trial = a0c + Yc.s
    bpost_c_trial = b0c + (nc-Yc.s)
    
    apost_c_hca = apost_c_trial + w*Ych
    bpost_c_hca = bpost_c_trial + w*(nch-Ych)
    
    apost_t = a0t + Yt.s 
    bpost_t = b0t + (nt-Yt.s)
    
    #calculate P(pt>pc|hybrid data) 
    P_ORRt_upper_Times_p_ORRc = function(y,ac,bc,at,bt){
      pbeta(y,at,bt,lower.tail=F)*dbeta(y,ac,bc)
    }
    
    pt_larger_pc[i] =integrate(P_ORRt_upper_Times_p_ORRc,
                                lower=0.0001, upper=0.9999,
                                ac = apost_c_hca, 
                                bc = bpost_c_hca,
                                at=apost_t,
                                bt=bpost_t)$value
  }
  #tau = quantile(pt_larger_pc, 1-alpha,type=1)
  tau = sort(pt_larger_pc)[round(nsim*(1-alpha))]
  
  return(tau)
}




