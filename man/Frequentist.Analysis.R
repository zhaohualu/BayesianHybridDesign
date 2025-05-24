#' Statistical Analysis Using Frequentist Method Based on Current Study Only
#'
#' This function performs the statistical analysis using frequentist method based on current study only.
#'
#' @param Yt Number of responses in experimental arm in current study
#' @param nt number of patients in experimental arm in current study
#' @param Yc Number of responses in control arm in current study
#' @param nc number of patients in control arm in current study
#' @param conflev Confidence Interval Level
#'
#' @return An object with values
#'  \itemize{
#'  \item pt Response rate of experimental arm
#'  \item pc Response rate of control arm
#'  \item delta Difference of response rate
#'  \item exactCI.c Clopper-Pearson CI for control arm
#'  \item exactCI.t Clopper-Pearson CI for experimental arm
#'  \item p.fisher Fishers exact test p-value. Always one-sided.
#'  }
#' @examples
#'
#' Frequentist.Analysis(Yt=18, nt=40,Yc=13, nc=40, conflev=0.8)
#'
#' @export
#'
Frequentist.Analysis= function(Yt=20, nt=40, Yc=12, nc=40, conflev=0.8){

  #Clopper-Pearson CI
  CI_c = exactci (r=Yc, n=nc, conflev=conflev)
  CI_t = exactci (r=Yt, n=nt, conflev=conflev)

  pt = Yt/nt
  pc = Yc/nc

  delta = pt - pc

  #Fisher's exact test
  p.fisher = fisher(Yc=Yc, nc=nc, Yt=Yt, nt=nt)$p.value

  return(list(pt=pt, pc=pc, delta=delta,
              exactCI.c = CI_c, exactCI.t = CI_t,
              p.fisher = p.fisher))
}

