#' Fishers Exact Test for Response Rate in Two arms
#' 
#' This function performs Fishers exact test customized for comparing response rate for two arms. 
#' 
#' @param r Number of reponses
#' @param n Total number of subjects
#' @param conflev Confidence level
#'  
#' @return Fishers test object
#' @examples
#' 
#' fisher(Yc=12, nc=40, Yt=19, nt=40)
#' 
#' @export
#' 
fisher = function(Yc=12, nc=40, Yt=19, nt=40, alternative = "greater"){
  #Yc: number of responders in control arm
  #nc: number of subjects in control arm
  #Yt: number of responders in experimental arm
  #nt: number of subjects in experimental arm
  #alternative: Default "greater"
  #Example: fisher(Yc=12, nc=40, Yt=19, nt=40)
  
  M <- as.table(rbind(c(Yt,Yc), c(nt-Yt, nc-Yc)))
  dimnames(M) <- list(Response = c("Response", "Non-resp"),
                      Group = c("Exp","control"))
  
  f = fisher.test(M, alternative = alternative)
  return(f)
}
