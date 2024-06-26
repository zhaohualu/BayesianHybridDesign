#' Support of Distributions
#'
#' Returns the support of a distribution.
#'
#' This function is exactly from the SAMprior R package version 1.1.1.
#' It is included to support the runSAM function in this package.
#'
#' takes x and transforms it according to the defined link function of
#' the mixture
#' @keywords internal
mixlink <- function(mix, x)
  attr(mix, "link")$link(x)

mixinvlink <- function(mix, x)
  attr(mix, "link")$invlink(x)

mixJinv_orig <- function(mix, x)
  attr(mix, "link")$Jinv_orig(x)

mixlJinv_orig <- function(mix, x)
  attr(mix, "link")$lJinv_orig(x)

mixlJinv_link <- function(mix, l)
  attr(mix, "link")$lJinv_link(l)

is.dlink <- function(x)
  inherits(x, "dlink")

is.dlink_identity <- function(x)
  is.dlink(x) & x$name == "identity"

is.mixidentity_link <- function(mix, l)
  is.dlink_identity(attr(mix, "link"))
#'
#' @param mix Mixture distribution.
#'
#' @keywords internal
support <- function(mix) UseMethod("support")
#' @export
support.default <- function(mix) "Unknown mixture"
#' @export
support.betaMix <- function(mix) mixlink(mix, c(0,1))
#' @export
support.gammaMix <- function(mix) mixlink(mix, c(0,Inf))
#' @export
support.normMix <- function(mix) mixlink(mix, c(-Inf,Inf))
