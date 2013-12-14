\name{yx_ssp.f}
\alias{BSplines}
\alias{CdN_DHxHt.f}
\alias{dN.f}
\alias{EYx_ssp.f}
\alias{Hx_root.f}
\alias{Int_CdN_DHx_dHt.f}
\alias{Int_E_VOL_AB_HmDm_HT_dHt.f}
\alias{qD.rout.f}
\alias{SK_EBLUP_LME.f}
\alias{SK_VOLab_EBLUP_LME.f}
\alias{TransKnots}
\alias{xy0_root.f}
\alias{xy0_SK_EBLUP_LME.f}
\alias{XZ_BSPLINE.f}
\alias{y2x_isp.f}
\alias{y2x_ssp.f}
\alias{yx_isp.f}
\alias{yx_ssp.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Internal funtions not usually called by users.
}
\description{
  Internal funtions not usually called by users.
}
\usage{
BSplines(knots = c(seq(0, 1, 0.1)), ord = 4, der = 0, x = c(seq(0, 1,
0.01)), ...)
CdN_DHxHt.f(Ht, Hx, qD, Hm, Dm, par.lme, ...)
dN.f(x, mw, sd, ...)
EYx_ssp.f(knt, coe, x, ...)
Hx_root.f(Hx, Dx, Hm, Dm, mHt, sHt, par.lme, ...)
Int_CdN_DHx_dHt.f(qD, Hx, Hm, Dm, mHt, sHt, par.lme, nGL = 51, ...)
Int_E_VOL_AB_HmDm_HT_dHt.f(Hm, Dm, A = NULL, B = NULL, iDH = "D",
mw_HtT, sd_HtT, par.lme, IA = F, nGL = 51, ...)
qD.rout.f(qD, alpha = 0.975, Hx, Hm, Dm, mHt, sHt, par.lme,
nGL = 51, ...)
SK_EBLUP_LME.f(xm, ym, xp, par.lme, ...)
SK_VOLab_EBLUP_LME.f(xm, ym, a = 0, b = 1, Ht, par.lme, IntPolOpt = T,
...)
TransKnots(knots = c(seq(0, 1, 0.1)), ord = 4, ...)
xy0_root.f(x, y0, SK, par.lme, ...)
xy0_SK_EBLUP_LME.f(xm, ym, y0, par.lme, ...)
XZ_BSPLINE.f(x, par.lme, ...)
y2x_isp.f(x, x.grd, y.grd, ...)
y2x_ssp.f(x, x.grd, y.grd, ...)
yx_isp.f(x, x.grd, y.grd, ...)
yx_ssp.f(x, x.grd, y.grd, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  %% we would need to describe the parameters of all alias-functions here...
  \item{knots}{
	Vector of knots
}
  \item{ord}{
	Scalar - order of spline
}
\item{der}{
  %%?
	der
}
\item{x}{
  %%?
	x
  }
    \item{Ht}{
	  %%     ??
	  Height
}
  \item{Hx}{
	%%     ~~Describe \code{Hx} here~~
	Numeric vector of stem heights (m) along which to return the expected diameter.
  }
    \item{Dx}{
	Scalar. Diameter for which to return height.
  }  
  \item{qD}{
%%     ~~Describe \code{qD} here~~
} 
  \item{Hm}{
	%%     ~~Describe \code{Hm} here~~
	Numeric vector of stem heights (m) along which diameter measurements were
	taken for calibration. Can be of length 1. Must be of same length as \code{Dm}. 
} 
  \item{Dm}{
	%%     ~~Describe \code{Dm} here~~
	Numeric vector of diameter measurements (cm)
	taken for calibration. Can be of length 1. Must be of same length as \code{Hm}.
  }  
  \item{par.lme}{
	%%     ~~Describe \code{par.lme} here~~
	List of taper model parameters obtained by \code{\link{TapeR_FIT_LME.f}}.
  }
  \item{mw}{
%%     ~~Describe \code{mw} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
} 
  \item{knt}{
%%     ~~Describe \code{knt} here~~
}
  \item{coe}{
%%     ~~Describe \code{coe} here~~
}
  \item{mHt}{
%%     ~~Describe \code{mHt} here~~
}
  \item{sHt}{
%%     ~~Describe \code{sHt} here~~
}
  \item{nGL}{
%%     ~~Describe \code{nGL} here~~
}
  \item{A}{
	Height from which to calculate
}
  \item{B}{
	Height to which to calculate
  }
    \item{iDH}{
%%     ~~Describe \code{iDH} here~~
}
  \item{mw_HtT}{
%%     ~~Describe \code{mw_HtT} here~~
}
  \item{sd_HtT}{
	%%     ~~Describe \code{sd_HtT} here~~
  }
  \item{IA}{
%%     ~~Describe \code{IA} here~~
}  
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{xm}{
%%     ~~Describe \code{xm} here~~
}
  \item{ym}{
%%     ~~Describe \code{ym} here~~
}
  \item{xp}{
%%     ~~Describe \code{xp} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{IntPolOpt}{
%%     ~~Describe \code{IntPolOpt} here~~
}
  \item{y0}{
%%     ~~Describe \code{y0} here~~
}
  \item{SK}{
%%     ~~Describe \code{SK} here~~
}
  \item{x.grd}{
%%     ~~Describe \code{x.grd} here~~
}
  \item{y.grd}{
%%     ~~Describe \code{y.grd} here~~
}
 \item{...}{Currently unused}
}
\details{
  Internal funtions not usually called by users.
  
  BSplines: Calculates BSplines%?

  CdN_DHxHt.f: %?

  dN.f:

  EYx_ssp.f:

  Hx_root.f:

  Int_CdN_DHx_dHt.f:

  Int_E_VOL_AB_HmDm_HT_dHt.f:
  
  qD.rout.f:
  
  SK_EBLUP_LME.f:
  
  SK_VOLab_EBLUP_LME.f:
  
  TransKnots:
  
  xy0_root.f:
  
  xy0_SK_EBLUP_LME.f:
  
  XZ_BSPLINE.f:
  
  y2x_isp.f:
  
  y2x_ssp.f:
  
  yx_isp.f:
  
  yx_ssp.f: 
}
\value{
  Depends on the specific function.
}
\author{
  %%  ~~who you are~~
  Edgar Kublin
}
% \note{
% %%  ~~further notes~~
% }
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{TapeR_FIT_LME.f}, \link{E_DHx_HmDm_HT.f}}
}
% \examples{
% ##
% }

