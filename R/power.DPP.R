#' Bayesian Hybrid Design
#'
#' This function calculates the power and design parameters Bayesian Hybrid Design using dynamic power prior approach.
#'
#' @param pt response rate for experimental arm in current study
#' @param nt number of patients in experimental arm in current study
#' @param pc Response rate for control arm in current study
#' @param nc number of patients in control arm in current study
#' @param pc.calib Required for calibration if tau is not provided. Response rate for control arm in current study for calibration. Usually, pc.calib = pch.
#' @param pch Response rate for control treatment in historical study
#' @param nche Equivalent number of patients borrowed from historical study
#' @param nch Total number of patients in historical control
#' @param alpha A scalar. One sided type I error rate. Required for calibration if tau is not provided.
#' @param tau Calibrated threshold for statistical significance. If tau is not provided, it will be calculated by calibration to type I error alpha.
#' @param a0c Hyperprior for control response rate beta(a0c, b0c)
#' @param b0c Hyperprior for control response rate beta(a0c, b0c)
#' @param a0t Hyperprior for experimental response rate beta(a0t, b0t)
#' @param b0t Hyperprior for experimental response rate beta(a0t, b0t)
#' @param delta_threshold Borrow when abs(pc_hat (current study) - pch) <= delta_threshold
#' @param method Method for dynamic borrowing, "Empirical Bayes", "Bayesian p", "Generalized BC", "JSD"
#' @param theta A parameter with a range of (0, 1), and applicable to method: "Generalized BC".
#' @param eta A parameter with a range of (0, infty), and applicable to method: "Bayesian p", "Generalized BC", "JSD". "Generalized BC" method requires two parameters theta and eta.
#' @param datamat A matrix with dimension nsim * 2 containing the pre-simulated data for the study treatment (1st column) and control  (1st column) groups, respectively. If not supplied, binomial random Monte Carlo samples will be generated in the function.
#' @param w0 prior power parameters w. If not specified (default), w_d is calculated by the specified method for dynamic borrowing.
#' @param nsim Number of replications to calculate power
#' @param seed=2000 seed for simulations
#'
#'
#'
#'
#' @return An object with values
#'  \itemize{
#'  \item power study power
#'  \item delta.bound Minimum detectable difference of response rate in the current study using hybrid design
#'  \item pc.PMD Posterior mean difference between hybrid control and concurrent control
#'  \item pc.sd  sd of PMD
#'  }
#' @examples
#'
#' o=power.DPP(pt=0.5, nt=40,pc=0.3,nc=40, pc.calib = 0.3, pch=0.3,nche=40,nch=180, alpha=0.1,a0c=0.001,b0c=0.001,a0t=0.001,b0t=0.001,delta_threshold=0.1)
#'
#' @export
#'
power.DPP = function(pt, nt, pc, nc, pc.calib,
                     pch,nche,nch,
                     alpha=0.1, tau=NULL,
                     a0c=0.001,b0c=0.001,a0t=0.001,b0t=0.001,
                     delta_threshold=0.1,
                     method="Empirical Bayes", theta=0.5, eta=1,
                     datamat = NULL, w0 = NULL,
                     nsim = 100000, seed=2000){
  #Threshold of significance tau for P(pt>pc|hybrid data) > tau
  if(is.null(tau)){
    tau = calibration (nt=nt, pc.calib=pc.calib, nc=nc, pch=pch,
                       nche=nche, nch=nch,
                       alpha = alpha, a0c=a0c, b0c=b0c, a0t=a0t, b0t=b0t,
                       delta_threshold=delta_threshold, method=method,
                       theta=theta, eta=eta, nsim = nsim, seed=seed)
  }

  #Dynamic borrowing parameter
  if(is.null(w0)){
    wt = borrow.wt(Yc=nc*pc, nc=nc, Ych=nch*pch, nch=nch, nche=nche,
                   a0c=a0c, b0c=b0c, delta_threshold=delta_threshold,
                   method=method,theta=theta,eta=eta)
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
  delta.bound = NA
  BoundaryIdx = which((phat_pt_larger_pc>=tau))
  if(length(BoundaryIdx)){
    delta.bound = (BoundaryIdx[1]-1)/nt-pc # - because SumDatat starts from 0
  }

  ######################
  # Power calculation
  ######################
  set.seed(seed)

  success = 0
  median_hca = median_c = rep(NA, nsim)
  mean_hca = mean_c = rep(NA, nsim)
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
                     method=method, theta=theta, eta=eta)
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
    # print("success")
    # print(success)
    # print("phat_pt_larger_pc")
    # print(phat_pt_larger_pc)
    # print("tau")
    # print(tau)

    phat_pt_larger_pc_all[i] = phat_pt_larger_pc

    success = success + (phat_pt_larger_pc>tau)

    median_hca[i] = qbeta(0.5, apost_c_hca, bpost_c_hca)
    median_c[i] = qbeta(0.5, apost_c_trial, bpost_c_trial)

    mean_hca[i] = apost_c_hca/(apost_c_hca + bpost_c_hca)
    mean_c[i] = apost_c_trial/(apost_c_trial + bpost_c_trial)

    wvec[i] = w

  }
  # print("success")
  # print(success)

  power = as.numeric(success) / nsim

  # print("power")
  # print(power)

  pc.PMD = mean(mean_hca-mean_c)
  pc.sd.PMD = sd(mean_hca-mean_c)

  return(list(power = power, tau = tau, alpha=alpha, pc.calib=pc.calib,
              nt=nt,nc=nc,nche=nche,pt=pt,pc=pc,nch=nch,pch=pch,
              method=method, theta=theta, eta=eta,
              delta_threshold=delta_threshold,
              pc.PMD=pc.PMD,pc.sd.PMD=pc.sd.PMD,
              delta.bound=delta.bound,
              phat_pt_larger_pc_all = phat_pt_larger_pc_all,
              mean_hca = mean_hca, mean_c = mean_c,
              simulated.data = datamat2,
              nsim=nsim, seed=seed,
              w=wvec))
}

