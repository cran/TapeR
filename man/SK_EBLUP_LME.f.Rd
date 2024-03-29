% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SK_EBLUP_LME.f.R
\name{SK_EBLUP_LME.f}
\alias{SK_EBLUP_LME.f}
\title{Evaluate fitted taper curve}
\usage{
SK_EBLUP_LME.f(xm, ym, xp, par.lme, Rfn = list(fn = "sig2"), ...)
}
\arguments{
\item{xm}{relative heights for which measurements are available}

\item{ym}{corresponding diameter measurements in height \code{xm}}

\item{xp}{relative heights for which predictions are required}

\item{par.lme}{Fitted model object, return of \code{\link{TapeR_FIT_LME.f}}}

\item{Rfn}{list with function name to provide estimated or assumed residual 
variances for the given measurements, optionally parameters for such functions}

\item{...}{not currently used}
}
\value{
a list holding nine elements:
\itemize{
  \item b_fix fixed effects parameter of taper model
  \item b_rnd random effects parameter given tree (posterior mean b_k)
  \item yp estimated diameter in height \code{xp}
  \item KOV_Mean variance-covariance matrix of expected value
  \item KOV_Pred variance-covariance matrix of prediction
  \item CI_Mean mean and limits of confidence interval
  \item MSE_Mean mean squared error of expected value
  \item MSE_Pred mean squared error of prediction
  \item CI_Pred mean and limits of prediction interval
}
}
\description{
This is the actual function to estimate diameters according to
the fitted mixed B-splines model.
}
\details{
This function is the actual working horse for prediction using the
fitted taper model. Based on the model \code{par.lme} and the measured
diameters \code{ym} and corresponding (relative) heights \code{xm} of a 
specific tree (there might be just one measurement), the random
effect parameters and subsequently diameters are estimated. Depending on the
parameter \code{Rfn}, the calibrated taper curve is forced through the
given diameter \code{ym} (\code{Rfn = list(fn="zero")}), or calibrated using
the complete residual variance-covariance information 
(\code{Rfn = list(fn="sig2")}, the default). 
Further assumptions are possible, see also \code{\link{resVar}} and
Kublin et al. (2013) p. 987 for more details.
}
\examples{
data("SK.par.lme")
TapeR:::SK_EBLUP_LME.f(1.3/27, 30, 1.3/27, SK.par.lme)
## using empirical best linear unbiased estimator: estimate != 30
TapeR:::SK_EBLUP_LME.f(1.3/27, 30, 1.3/27, SK.par.lme, Rfn=list(fn="sig2"))$yp
## interpolate / force through given diameter: estimate  == 30
TapeR:::SK_EBLUP_LME.f(1.3/27, 30, 1.3/27, SK.par.lme, Rfn=list(fn="zero"))$yp
TapeR:::SK_EBLUP_LME.f(1.3/27, 30, c(1.3, 5)/27, SK.par.lme)
par.lme <- SK.par.lme
h <- 12 # tree height
xm <- c(1.3, 3) / h # relative measuring height
ym <- c(8, 7.5) # measured diameter
xp <- c(0.5, 1) / h # relative prediction height
TapeR:::SK_EBLUP_LME.f(xm, ym, xp, SK.par.lme)
}
\seealso{
\code{\link{E_DHx_HmDm_HT.f}}, \code{\link{E_VOL_AB_HmDm_HT.f}}, 
\code{\link{resVar}}
}
\author{
Edgar Kublin
}
