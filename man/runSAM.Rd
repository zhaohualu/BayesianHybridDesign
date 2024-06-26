% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runSAM.R
\name{runSAM}
\alias{runSAM}
\alias{runSAM.betaMix}
\alias{runSAM.normMix}
\title{Generating Operating Characteristics of SAM Priors}
\usage{
runSAM(
  if.prior,
  theta.h,
  method.w,
  prior.odds,
  nf.prior,
  delta,
  n,
  n.t,
  decision,
  ntrial,
  if.MAP,
  weight,
  theta,
  theta.t,
  datamat = NULL,
  ...
)

\method{runSAM}{betaMix}(
  if.prior,
  theta.h,
  method.w,
  prior.odds,
  nf.prior,
  delta,
  n,
  n.t,
  decision,
  ntrial,
  if.MAP,
  weight,
  theta,
  theta.t,
  datamat = NULL,
  ...
)

\method{runSAM}{normMix}(
  if.prior,
  theta.h,
  method.w,
  prior.odds,
  nf.prior,
  delta,
  n,
  n.t,
  decision,
  ntrial,
  if.MAP,
  weight,
  theta,
  theta.t,
  datamat = NULL,
  ...,
  sigma
)
}
\arguments{
\item{if.prior}{Informative prior constructed from historical data,
represented (approximately) as a mixture of conjugate distributions.}

\item{theta.h}{Estimate of the treatment effect based on historical data.
If missing, the default value is set to be the posterior mean estimate from
\code{if.prior}.}

\item{method.w}{Methods used to determine the mixture weight for SAM priors.
The default method is LRT (Likelihood Ratio Test), the alternative option can
be PPR (Posterior Probability Ratio). See \code{\link{SAM_weight}} for more
details.}

\item{prior.odds}{The prior probability of \eqn{H_0} being true compared to
the prior probability of \eqn{H_1} being true using PPR method. The default
value is 1. See \code{\link{SAM_weight}} for more details.}

\item{nf.prior}{Non-informative prior used for constructing the SAM prior
and robust MAP prior.}

\item{delta}{Clinically significant difference used for the SAM prior.}

\item{n}{Sample size for the control arm.}

\item{n.t}{Sample size for the treatment arm.}

\item{decision}{Decision rule to compare the treatment with the control;
see \code{\link{decision2S}}.}

\item{ntrial}{Number of trials simulated.}

\item{if.MAP}{Whether to simulate the operating characteristics of the
robust MAP prior for comparison, the default value is \code{FALSE}.}

\item{weight}{Weight assigned to the informative prior component
(\eqn{0 \leq} \code{weight} \eqn{\leq 1}) for the robust MAP prior,
the default value is 0.5.}

\item{theta}{A vector of the response rate (binary endpoints) or mean
(continuous endpoints) for the control arm.}

\item{theta.t}{A vector of the response rate (binary endpoints) or mean
(continuous endpoints) for the treatment arm.}

\item{datamat}{A ntrial * 2 matrix of simulated binomial samples,
first column is study treament group treatment, second column is study control group}

\item{...}{Additional parameters for continuous endpoints.}

\item{sigma}{Variance to simulate the continuous endpoint under normality
assumption.}
}
\value{
Returns a list
result - dataframe that contains the relative bias, relative MSE,
type I error, and power for both SAM priors, as well as robust MAP priors.
Additionally, the mixture weight of the SAM prior is also displayed.

simulated.data - a matrix of two columns, first treatment, second control
post_theta_t_list - list of all replication, posterior distribution of treatment group
post_theta_c_list - list of all replication, posterior distribution of control group, non informative prior
post_theta_c_SAM_list - list of all replication, posterior distribution of control group, SAM prior
post_theta_c_MAP_list - list of all replication, posterior distribution of control group, MAP prior
}
\description{
This function is modified based on the get_OC function from the SAMprior R package version 1.1.1.
}
\details{
This function is modified based on the \code{get_OC} from the SAMprior R package version 1.1.1.
The modification is made in order that the function can take simulated data (the datamat argument) from outside the function rather than simulating data within the funciton
We used the same data for all methods for comparison in order to get a better comparison.

The \code{runSAM} function is designed to generate the operating
characteristics of SAM priors (\emph{Yang, et al., 2023}), including the
relative bias, relative mean squared error, and type I error and power
under a two-arm comparative trial design. As an option, the operating
characteristic of robust MAP priors (\emph{Schmidli, et al., 2014})
can also be generated for comparison.

The \code{runSAM} function is designed to generate the operating
characteristics of SAM priors, including the relative bias, relative
mean squared error, and type I error, and power under a two-arm
comparative trial design. As an option, the operating characteristics of
robust MAP priors (\emph{Schmidli, et al., 2014}) can also be generated for
comparison.

The relative bias is defined as the difference between the bias of a method
and the bias of using a non-informative prior. The relative mean squared
error is the difference between the mean squared error (MSE) of a method and
the MES of using a non-informative prior.

To evaluate type I error and power, the determination of whether the
treatment is superior to the control is calculated based on function
\code{\link{decision2S}}.
}
\section{Methods (by class)}{
\itemize{
\item \code{runSAM(betaMix)}: The function is designed to generate the operating
characteristics of SAM priors for binary endpoints.

\item \code{runSAM(normMix)}: The function is designed to generate the operating
characteristics of SAM priors for continuous endpoints.

}}
\examples{
set.seed(123)
## Example of a binary endpoint
## Consider a randomized comparative trial designed to borrow information
## from historical data on the control. We assumed a non-informative prior
## beta(1, 1) and an informative prior beta(30, 50) after incorporating
## the historical data. The treatment is regarded as superior to the control
## if Pr(RR.t > RR.c | data) > 0.95, where RR.t and RR.c are response rates
##  of the treatment and control, respectively. The operating characteristics
##  were assessed under the scenarios of (RR.c, RR.t) = (0.3, 0.36) and (0.3, 0.56).
## OC <- runSAM(## Informative prior constructed based on historical data
##              if.prior = mixbeta(c(1, 30, 50)),
##              ## Non-informative prior used for constructing the SAM prior
##              nf.prior = mixbeta(c(1,1,1)),
##              delta    = 0.2,  ## Clinically significant difference
##              n = 35,          ## Sample size for the control arm
##              n.t = 70,        ## Sample size for the treatment arm
##              ## Decision rule to compare the whether treatment is superior
##              ## than the control
##              decision = decision2S(0.95, 0, lower.tail=FALSE),
##              ntrial   = 1000,  ## Number of trials simulated
##              ## Weight assigned to the informative component for MAP prior
##              weight = 0.5,
##              ## A vector of response rate for the control arm
##              theta    = c(0.3, 0.36),
##              ## A vector of response rate for the treatment arm
##              theta.t  = c(0.3, 0.56))
## OC

## Example of continuous endpoint
## Consider a randomized comparative trial designed to borrow information
## from historical data on the control. We assumed a non-informative prior
## N(0, 1e4) and an informative prior N(0.5, 2) after incorporating
## the historical data. The treatment is regarded as superior to the control
## if Pr(mean.t > mean.c | data) > 0.95, where mean.t and mean.c are mean
##  of the treatment and control, respectively. The operating characteristics
##  were assessed under the scenarios of (mean.c, mean.t) = (0.1, 0.1) and
## (0.5, 1.0).
sigma      <- 2
prior.mean <- 0.5
prior.se   <- sigma/sqrt(100)
## OC <- runSAM(## Informative prior constructed based on historical data
##              if.prior = mixnorm(c(1, prior.mean, prior.se)),
##              ## Non-informative prior used for constructing the SAM prior
##              nf.prior = mixnorm(c(1, 0, 1e4)),
##              delta    = 0.2 * sigma,  ## Clinically significant difference
##              n = 100,                 ## Sample size for the control arm
##              n.t = 200,               ## Sample size for the treatment arm
##              ## Decision rule to compare the whether treatment is superior
##              ## than the control
##              decision = decision2S(0.95, 0, lower.tail=FALSE),
##              ntrial   = 1000,  ## Number of trials simulated
##              ## A vector of mean for the control arm
##              theta    = c(0.1, 0.5),
##              ## A vector of mean for the treatment arm
##              theta.t  = c(0.1, 1.0),
##              sigma = sigma)
## OC

}
\references{
Yang P, Zhao Y, Nie L, Vallejo J, Yuan Y.
SAM: Self-adapting mixture prior to dynamically borrow information from
historical data in clinical trials. \emph{Biometrics} 2023; 00, 1–12.
https://doi.org/10.1111/biom.13927

Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B.
Robust meta-analytic-predictive priors in clinical trials with historical control information.
\emph{Biometrics} 2014; 70(4):1023-1032.
}
