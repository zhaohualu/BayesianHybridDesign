#' Rejection boundary by Fisher's exact test
#'
#' This function calculates the minimum detectable difference in response rate when using Fisher's exact test.
#'
#' @param pc response rate for control arm
#' @param nc number of patients in control arm
#' @param nt number of patients in experimental arm
#' @param alpha p-value threshold for significance. Alpha must be one-sided. Default 0.1.
#'
#' @return An object with values
#'  \itemize{
#'  \item M A 2-by-2 contingency table for response and treatment group
#'  \item p p-value when the response rates just across the boundary
#'  \item rc number of responders for control arm, determined by the response rate and sample size
#'  \item nc sample size in control arm
#'  \item rt response rate for experimental arm just acros the boundary
#'  \item nt sample size in experimental arm
#'  \item delta Minimum detectable difference of response rate between two groups
#'  }
#' @examples
#'
#' fisher.bound(pc=0.3,nc=40,nt=40,alpha=0.1)
#'
#' @export
#'
fisher.bound = function(pc,nc,nt,alpha=0.1){
  rc = round(nc*pc)

  for(rt in 1:nt){
    M <- as.table(cbind(c(rt,nt-rt),c(rc,nc-rc)))
    dimnames(M) <- list(Response = c("Response", "Non-resp"),
                        Group = c("Exp","control"))

    o = fisher.test(M,alternative = "greater")

    p = o$p.value
    if(p <= alpha){
      break
    }
  }

  deltaboundary = rt/nt-rc/nc

  return(list(M=M,
              p = p,
              rc = rc,
              nc = nc,
              rt = rt,
              nt = nt,
              delta = deltaboundary))
}
