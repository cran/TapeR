\name{TapeR-package}
\alias{TapeR-package}
\alias{TapeR}
\docType{package}
\title{
  Flexible tree taper curves based on Semiparametric Mixed Models.
}
\description{
  This package implements functions for fitting taper curves (a semiparametric linear
  mixed effects taper model) to diameter measurements along
  stems. Further functions are provided to estimate the
  variance/confidence intervals
  around the predicted curves, to calculate timber volume (also by
  sections) and marginal (e.g., upper) diameters. For cases where tree
  heights are not measured, methods for estimating
  additional variance in volume predictions resulting from uncertainties in
  tree height models (tariffs) are provided. The
  example data include the taper curve parameters for Norway spruce used
  in the 3rd German NFI fitted to 380 trees and a subset of section-wise
  diameter measurements of these trees.
}
\details{
\tabular{ll}{
Package: \tab TapeR\cr
Type: \tab Package\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Fits taper models using diameter measurements along the stem. Uses fitted
models and arbitrary numbers of diameter measurements to estimate
diameter at any position along the stem. Estimates timber volume from the
taper curve. Provides variances for all estimates.
}
\author{
  Edgar Kublin
  
  Maintainer: Johannes Breidenbach <job@skogoglandskap.no>
}

\references{
Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem taper
and volume prediction method based on mixed-effects B-spline regression,
Eur J For Res, 132:983-997. 
}

\keyword{ package }
\seealso{
  \code{\link{TapeR_FIT_LME.f}}, \code{\link{E_DHx_HmDm_HT.f}},
  \code{\link{DxHx.df}}, \code{\link{SK.par.lme}}, \code{\link{HT.par}}
}
