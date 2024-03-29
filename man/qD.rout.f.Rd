% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qD.rout.f.R
\name{qD.rout.f}
\alias{qD.rout.f}
\title{qD.rout.f}
\usage{
qD.rout.f(qD, alpha = 0.975, Hx, Hm, Dm, mHt, sHt, par.lme, Rfn, nGL = 51, ...)
}
\arguments{
\item{qD}{vector of quantiles, finally passed to \code{pnorm}}

\item{alpha}{quantile for which root is sought}

\item{Hx}{Numeric vector of stem heights (m) along which to return the 
expected diameter}

\item{Hm}{Numeric vector of stem heights (m) along which diameter measurements 
were taken for calibration. Can be of length 1. Must be of same length as \code{Dm}}

\item{Dm}{Numeric vector of diameter measurements (cm) taken for calibration.
Can be of length 1. Must be of same length as \code{Hm}}

\item{mHt}{Scalar. Tree height (m)}

\item{sHt}{Scalar. Standard deviation of stem height. Can be 0 if height was 
measured without error}

\item{par.lme}{List of taper model parameters obtained by 
\code{\link{TapeR_FIT_LME.f}}.}

\item{Rfn}{list with function name to provide estimated or assumed residual 
variances for the given measurements, optionally parameters for such functions}

\item{nGL}{Numeric scalar. Number of support points for numerical integration.}

\item{...}{not currently used}
}
\value{
\code{qD} for given \code{alpha} with respect to \code{Int_CdN_DHx_dHt.f}
}
\description{
Internal function not usually called by users
}
\author{
Edgar Kublin
}
