#' Plot of DPP Object produced by DPP.analysis() function.
#' 
#' Plot the posterior distributions of OR in control arm and experimental arm. DPP object produced by DPP.analysis() function. 
#'
#' This function Perform Analysis for a Study Using Bayesian Hybrid Design using Dynamic Power Prior Framework.
#'
#' @param  DPP DPP object produced by DPP.analysis() function. 
#' 
#' @examples
#'
#' o = DPP.analysis(Yt=39, nt=60, Yc=13, nc=30, Ych=90, nch=200, nche = 30, 
#' a0c= 0.001, b0c= 0.001, a0t= 0.001, b0t= 0.001,
#' delta_threshold = 0.1, method = "Empirical Bayes", 
#' theta = 0.5, eta = 1)
#' 
#' plot.DPP(DPP=o)
#' 
#' @export
#' 
plot.DPP = function(DPP){
  ### plot the posterior density curves
  x = seq(0, 1, length=100)
  yc = dbeta(x, DPP$apost_c_trial, DPP$bpost_c_trial)
  yhc = dbeta(x, DPP$apost_c_hca, DPP$bpost_c_hca)
  yt = dbeta(x, DPP$apost_t, DPP$bpost_t)
  
  # Plot the first density curve
  plot(x, yc, col = "blue", lwd = 2, type="l",
       ylim = c(0, max(yc, yhc, yt)), 
       main = paste("Posterior ORR: P(ORR_t > ORR_c): ", round(DPP$phat_pt_larger_pc, 3)), xlab = "ORR", ylab = "Density")
  
  # Add the second density curve
  lines(x, yhc, col = "red", lwd = 2)
  lines(x, yt, col = "green", lwd = 2)
  
  # Fill under first curve
  polygon(c(x, rev(x)), c(yc, rep(0, length(yc))), 
          col = rgb(0, 0, 1, 0.3), border = NA)
  
  # Fill under second curve
  polygon(c(x, rev(x)), c(yhc, rep(0, length(yhc))), 
          col = rgb(1, 0, 0, 0.3), border = NA)
  
  
  # Fill under second curve
  polygon(c(x, rev(x)), c(yhc, rep(0, length(yhc))), 
          col = rgb(1, 0, 0, 0.3), border = NA)
  
  # Fill under experimental curve
  polygon(c(x, rev(x)), c(yt, rep(0, length(yt))), 
          col = rgb(0, 1, 0, 0.3), border = NA)
  
  # Add legend
  legend("topleft", legend = c("Concurrent Control", "Hybrid Control", "Experimental Arm"), 
         fill = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.3), rgb(0, 1, 0, 0.3)), bty = "n", cex = 0.7)
}

