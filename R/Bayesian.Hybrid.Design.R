#' Bayesian Hybrid Design
#' 
#' This function calculates the power and design parameters Bayesian Hybrid Design using dynamic power prior approach.
#' 
#' @param pt response rate for experimental arm in current study
#' @param nt number of patients in experimental arm in current study    
#' @param nc number of patients in control arm in current study  
#' @param nche Equivalent number of patients borrowed from historical study
#' @param nch Total number of patients in historical control
#' @param pc Response rate for control arm in current study
#' @param pch Response rate for control treatment in historical study
#' @param sig Significance boundary defined as P(pt_hat > pc_hat|hybrid) > sig
#' @param a0c Hyperprior for control response rate beta(a0c, b0c)
#' @param b0c Hyperprior for control response rate beta(a0c, b0c)
#' @param a0t Hyperprior for experimental response rate beta(a0t, b0t)
#' @param b0t Hyperprior for experimental response rate beta(a0t, b0t)
#' @param delta_threshold Borrow when abs(pc_hat (current study) - pch) <= delta_threshold 
#' @param method Method for dynamic borrowing, "Empirical Bayes" or "Heterogeneity".
#' @param datamat A matrix with dimension nsim * 2 containing the pre-simulated data for the study treatment (1st column) and control  (1st column) groups, respectively.
#' @param nsim Number of replications to calculate power
#' @param seed=2000 seed for simulations
#' @param method Method for dynamic borrowing, "Empirical Bayes", "Bayesian p", "Generalized BC", "JSD"
#' 
#' @return An object with values
#'  \itemize{
#'  \item power study power
#'  \item delta.bound Minimum detectable difference of response rate in the current study using hybrid design
#'  \item pc.bias bias estimate of posterior median ORR for control arm using hybrid design
#'  \item pc.mse MSE of posterior median ORR for control arm using hybrid design
#'  }
#' @examples
#' 
#' Bayesian.Hybrid.Design(pt=0.512, nt=40,pc=0.312,nc=40,pch=0.312,nche=40,nch=234, sig=0.90,a0c=0.001,b0c=0.001,a0t=0.001,b0t=0.001,delta_threshold=0.1)
#' 
#' @export
#' 
Bayesian.Hybrid.Design = function(pt=0.512, nt=40,pc=0.312,nc=40,pch=0.312,
                                  nche=40,nch=234, sig=0.90,
                                  a0c=0.001,b0c=0.001,a0t=0.001,b0t=0.001,
                                  delta_threshold=0.1, 
                                  method="Empirical Bayes", theta=0.5,
                                  datamat = NULL, w0 = NULL,
                                  nsim = 10000, seed=2000){
  
  #Dynamic borrowing parameter
  if(is.null(w0)){
    wt = borrow.wt(Yc=nc*pc, nc=nc, Ych=nch*pch, nch=nch, nche=nche, 
                   a0c=a0c, b0c=b0c, delta_threshold=delta_threshold, method=method)
    w = wt$w
  } else {
    w=w0
  }
  #All scenarios for number of responders in exp arm
  Yt = 0:nt # start from 0
  Yc = nc*pc #number of responders in control arm
  
  #number of responders in historical control arm
  Ych=nch*pch
  
  #Poster distributions
  apost_c_trial = a0c + Yc
  bpost_c_trial = b0c + (nc-Yc)
  
  apost_c_hca = apost_c_trial + w*Ych
  bpost_c_hca = bpost_c_trial + w*(nch-Ych)
  
  apost_t = a0t + Yt 
  bpost_t = b0t + (nt-Yt)
  
  #Boundary number of responders (BoundaryIdx-1) for significance
  phat_pt_larger_pc = rep(NA, nt+1)
  
  for(i in 1:length(Yt)){
    P_ORRt_upper_Times_p_ORRc = function(y,ac,bc,at,bt){
      pbeta(y,at,bt,lower.tail=F)*dbeta(y,ac,bc)
    }
    phat_pt_larger_pc[i]=integrate(P_ORRt_upper_Times_p_ORRc,
                                   lower=0.0001, upper=0.9999,
                                   ac = apost_c_hca, bc = bpost_c_hca,
                                   at=apost_t[i],bt=bpost_t[i])$value
  }
  
  #min detectable response difference: delta.bound
  BoundaryIdx = which((phat_pt_larger_pc>=sig))
  if(length(BoundaryIdx)){
    delta.bound = (BoundaryIdx[1]-1)/nt-pc # - because SumDatat starts from 0
  }
  
  ######################
  # Power calculation
  ######################
  set.seed(seed)
  
  success = 0
  median_hca = median_c = rep(NA, nsim)
  phat_pt_larger_pc_all = rep(NA, nsim)
  
  datamat2 = matrix(NA,nsim,2)
  
  wvec = numeric(nsim)
  
  #for each simulated trial, determine whether significant  
  for(i in 1:nsim){
    if(is.null(datamat)){
      Yt.s = rbinom(1,size=nt,prob=pt) #simulated Yt
      Yc.s = rbinom(1,size=nc,prob=pc) #simulated Yc
    } else{
      Yt.s =datamat[i,1]
      Yc.s =datamat[i,2]
    }
    
    datamat2[i,] = c(Yt.s,Yc.s)
    
    #Dynamic borrowing parameter
    if(is.null(w0)){
      wt = borrow.wt(Yc=Yc.s, nc=nc, Ych=nch*pch, nch=nch, nche=nche,  
                     a0c=a0c, b0c=b0c, delta_threshold=delta_threshold, 
                     method=method, theta=theta)
      w = wt$w
    }else{
      w = w0
    }
    #Posterior distributions
    apost_c_trial = a0c + Yc.s
    bpost_c_trial = b0c + (nc-Yc.s)
    
    apost_c_hca = apost_c_trial + w*Ych
    bpost_c_hca = bpost_c_trial + w*(nch-Ych)
    
    apost_t = a0t + Yt.s 
    bpost_t = b0t + (nt-Yt.s)
    
    #calculate the probability of 
    P_ORRt_upper_Times_p_ORRc = function(y,ac,bc,at,bt){
      pbeta(y,at,bt,lower.tail=F)*dbeta(y,ac,bc)
    }
    
    phat_pt_larger_pc=integrate(P_ORRt_upper_Times_p_ORRc,
                                lower=0.0001, upper=0.9999,
                                ac = apost_c_hca, 
                                bc = bpost_c_hca,
                                at=apost_t,
                                bt=bpost_t)$value
    
    phat_pt_larger_pc_all[i] = phat_pt_larger_pc
    
    success = success + (phat_pt_larger_pc>sig)
    
    median_hca[i] = qbeta(0.5, apost_c_hca, bpost_c_hca)
    median_c[i] = qbeta(0.5, apost_c_trial, bpost_c_trial)
    
    wvec[i] = w
    
  }
  power = success / nsim
  pc.bias = mean(median_hca-median_c)
  pc.mse = mean((median_hca-median_c)^2)
  
  return(list(power = power, 
              delta.bound=delta.bound, 
              pc.bias=pc.bias, 
              pc.mse=pc.mse, 
              phat_pt_larger_pc_all = phat_pt_larger_pc_all,
              median_hca = median_hca,
              median_c = median_c,
              simulated.data = datamat2,
              nsim=nsim, seed=seed,
              w=wvec))
}

