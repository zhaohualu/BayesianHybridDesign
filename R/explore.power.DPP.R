#'  Power Calculation of Several power.DPP Function Parameters in Bayesian Hybrid Design
#'
#' This function calculates the power and design parameters Bayesian Hybrid Design using
#' parallel computing for efficiency. It is useful when there are multiple values of X
#' parameters to explore. If only one setting to calculate, please use the power.DPP function.
#' What is exploring, what parameters allow vectors for exploration.
#'
#' @param method A vector of methods for dynamic borrowing for exploration, supports "Empirical Bayes", "Bayesian p", "Generalized BC", "JSD"
#' @param pc A vector of response rate for study control arm for exploration
#' @param nc A vector of sample size for study control arm for exploration
#' @param pc.calib A vector of response rate of study control arm for exploration. Each of these response rates is used for calibration such that the p-value equals to alpha (specified type I error) if pc = pc.calib
#' @param q A vector of ratio of historical control arm and the study control arm for exploration.
#'
#' @param delta A scalar of difference between pt and pc to determine pt in each setting.
#' @param r A scalar of  ratio of study experimental arm and study control arm  for exploration
#' @param pch A scalar of the response rate for historical control arm
#' @param nch A scalar of the total number of patients in historical control
#' @param alpha A scalar. One sided type I error rate.
#' @param a0c A scalar. Hyperprior for control response rate beta(a0c, b0c). Default 0.001.
#' @param b0c A scalar. Hyperprior for control response rate beta(a0c, b0c). Default 0.001.
#' @param a0t A scalar. Hyperprior for experimental response rate beta(a0t, b0t). Default 0.001.
#' @param b0t A scalar. Hyperprior for experimental response rate beta(a0t, b0t). Default 0.001.
#' @param delta_threshold Borrow when abs(pc_hat (current study) - pch) <= delta_threshold#'. Default 0.1.
#' @param theta A scalar parameter with a range of (0, 1), and applicable to method: "Generalized BC". Default 0.5.
#' @param eta A scalar parameter with a range of (0, infty), and applicable to method: "Bayesian p", "Generalized BC", "JSD". "Generalized BC" method requires two parameters theta and eta. Default 1.
#' @param nsim A scalar. Number of replications to calculate power. Default 100,000.
#' @param seed A scalar. seed for simulations. Default 2000.
#' @param ncore description
#'
#'
#'
#'
#' @return A dataframe containing the setting parameters and
#'  \itemize{
#'  \item type I error
#'  \item power
#'  \item delta boundary
#'  \item mean PMD
#'  \item sd(PMD)
#'  }
#' @examples
#'
#' Res1 = explore.power.DPP(method=c("Empirical Bayes"),
#' pc=c(0.17,0.27,0.37),
#' nc=23:33,
#' pc.calib = 0.27,
#' q = c(1,1.5,2),
#' pch=0.27, nch=500,
#' alpha=0.1,
#' a0c=0.001, b0c=0.001, a0t=0.001, b0t=0.001,
#' delta_threshold=0.1,
#' theta=0.5, eta=1,
#' nsim = 1000, seed=2000,
#' ncore=16)
#'
#' @export
#'
#' @import parallel
#' @import doParallel
#' @import foreach
explore.power.DPP= function(method,
                            pc,
                            nc,
                            pc.calib,
                            q,
                            delta,
                            r,
                            pch,
                            nch,
                            alpha = 0.1,
                            a0c=0.001,
                            b0c=0.001,
                            a0t=0.001,
                            b0t=0.001,
                            delta_threshold=0.1,
                            theta=0.5,
                            eta=1,
                            nsim = 100000,
                            seed=2000,
                            ncore=NULL){

  # Setting up the parallel computing backend
  if(is.null(ncore)){
    ncore = detectCores()-1
  }
  cl <- makeCluster(ncore)
  registerDoParallel(cl)

  # Generate the setting matrix, each row is one setting
  Settings = expand.grid(method = method, pc = pc, nc = nc, pc.calib = pc.calib, q = q)

  Settings$nt = Settings$nc * r
  Settings$nche = round(Settings$nc * Settings$q)

  # Loop over all settgins with parallel computing
  NS = nrow(Settings)
  Result = foreach(i = 1:NS,.combine=rbind) %dopar%{

    # Getting setting parameters for each setting
    pci = Settings$pc[i]
    nci = Settings$nc[i]
    methodi = Settings$method[i]
    nti = Settings$nt[i]
    nchei = Settings$nche[i]
    pc.calibi = Settings$pc.calib[i]
    qi =  Settings$q[i]

    pti = pci + delta

    # Get the calibrated threshold tau
    # Get power, delta boundary, mean PMD and sd(PMD)
    o = power.DPP (pt=pti, nt=nti, pc=pci, nc=nci, pc.calib = pc.calibi,
                   pch=pch, nche=nchei,nch=nch,
                   alpha=alpha, tau=NULL,
                   a0c=a0c, b0c=b0c, a0t=a0t, b0t=b0t,
                   delta_threshold=delta_threshold, method=methodi,
                   theta=theta, eta=eta, nsim = nsim, seed=seed)
    tau = o$tau
    pow = o$power
    delta.bound=o$delta.bound
    mean.PMD=o$pc.PMD
    sd.PMD=o$pc.sd.PMD

    mean_pc_hca = mean(o$mean_hca)
    mean_pc_c = mean(o$mean_c)

    # Based on the calibrated tau, calculate the type I error given pt = pc in the current setting
    o = power.DPP (pt=pci, nt=nti, pc=pci, nc=nci, pc.calib = pc.calibi,
                   pch=pch,nche=nchei,nch=nch,
                   alpha=alpha, tau=tau,
                   a0c=a0c, b0c=b0c, a0t=a0t, b0t=b0t,
                   delta_threshold=delta_threshold, method=methodi,
                   theta=theta, eta=eta, nsim = nsim, seed=seed)
    type1err = o$power

    # Return the setting parameters and type I error, power, delta boundary, mean PMD and sd(PMD)
    ans = data.frame(method = methodi,
                     nt = nti,
                     nc = nci,
                     nche = nchei,
                     pt = pti,
                     pc = pci,
                     pc.calib = pc.calibi,
                     pch = pch,
                     tau = tau,
                     type1err = type1err,
                     power = pow,
                     delta.boundary = delta.bound,
                     mean.PMD = mean.PMD,
                     sd.PMD = sd.PMD,
                     q = qi)

    return(ans)
  }

  # Clearn up the parallel computing environment
  stopCluster(cl)

  return(Result)



  # Result2$delta = deltai
  # Result2$method = methodi
}



