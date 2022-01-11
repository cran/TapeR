#' @title Fits a taper curve model to the specified diameter-height data
#' @description Fits a taper curve model with random effects on tree-level based
#' on B-Splines to the specified diameter-height data. Number and position of
#' nodes and order of B-Splines can be specified.
#' @param Id Vector of tree identifiers of same length as diameter and height
#' measurements.
#' @param x Numeric vector of height measurements (explanatory variables) along
#' the stem relative to the tree height.
#' @param y Numeric vector of diameter measurements (response) along the stem
#' (in centimeters).
#' @param knt_x Numeric vector of relative knot positions for fixed effects.
#' @param ord_x Numeric scalar. Order of fixed effects Spline (4=cubic).
#' @param knt_z Numeric vector of relative knot positions for random effects.
#' @param ord_z Numeric scalar. Order of random effects Spline (4=cubic).
#' @param IdKOVb Character string. Type of covariance matrix used by
#' \code{lme}. Only "pdSymm" makes sense. Rather reduce number of knots if
#' function does not converge.
#' @param control a list of control values for the estimation algorithm to
#' replace the default values returned by the function
#' \code{\link[nlme]{lmeControl}}. Defaults to an empty list.
#' @param ... not currently used
#' @details If too few trees are given, the linear mixed model (lme) will not
#' converge. See examples for a suggestion of node positions.
#'
#' The variance parameters \code{theta} are stored in the natural parametrization
#' (Pinheiro and Bates (2004), p. 93). This means log for variances and logit for
#' covariances. \code{theta} is the vectorized triangle of the random effects
#' covariance matrix + the residual variance (lSigma). Given there are 2 inner
#' knots for random effects, the structure will be
#' c(sig^2_b1, sig_b1 sig_b2, sig_b1 sig_b3, sig_b1 sig_b4, sig^2_b2,...,sig^2_b4, lSigma)
#' @return List of model properties
#' \itemize{
#'   \item{fit.lme}{Summary of the fitted lme model.}
#'   \item{par.lme}{List of model parameters (e.g., coefficients and
#'   variance-covariance matrices) needed for volume estimation and other
#'   functions in this package.}
#'   Components of the \code{par.lme} list
#'   \itemize{
#'     \item{knt_x}{Relative positions of the fixed effects Spline knots along the stem.}
#'     \item{pad_knt_x}{Padded version of knt_x, as used to define B-Spline design matrix.}
#'     \item{ord_x}{Order of the spline.}
#'     \item{knt_z}{Relative positions of the random effects Spline knots along
#'     the stem.}
#'     \item{pad_knt_z}{Padded version of knt_z, as used to define B-Spline design matrix.}
#'     \item{ord_z}{Order of the spline.}
#'     \item{b_fix}{Fixed-effects spline coefficients.}
#'     \item{KOVb_fix}{Covariance of fixed-effects.}
#'     \item{sig2_eps}{Residual variance.}
#'     \item{dfRes}{Residual degrees of freedom.}
#'     \item{KOVb_rnd}{Covariance of random effects.}
#'     \item{theta}{Variance parameters in natural parametrization. See Details. }
#'     \item{KOV_theta}{Approximate asymptotic covariance matrix of variance parameters.}
#'   }
#' }
#' @author Edgar Kublin
#' @references Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem
#' taper and volume prediction method based on mixed-effects B-spline
#' regression, Eur J For Res, 132:983-997.
#' @seealso \code{\link{E_DHx_HmDm_HT.f}}, \code{\link{E_DHx_HmDm_HT_CIdHt.f}},
#' \code{\link{E_HDx_HmDm_HT.f}}, \code{\link{E_VOL_AB_HmDm_HT.f}}
#' @import nlme
#' @importFrom stats vcov cor cov2cor anova
#' @export
#'
#' @examples
#' # load example data
#' data(DxHx.df)
#'
#' # prepare the data (could be defined in the function directly)
#' Id = DxHx.df[,"Id"]
#' x = DxHx.df[,"Hx"]/DxHx.df[,"Ht"]#calculate relative heights
#' y = DxHx.df[,"Dx"]
#'
#' # define the relative knot positions and order of splines
#' knt_x = c(0.0, 0.1, 0.75, 1.0);	ord_x = 4 # B-Spline knots: fix effects; order (cubic = 4)
#' knt_z = c(0.0, 0.1       ,1.0); ord_z = 4 # B-Spline knots: rnd effects
#'
#' # fit the model
#' taper.model <- TapeR_FIT_LME.f(Id, x, y, knt_x, ord_x, knt_z, ord_z,
#'                                IdKOVb = "pdSymm")
#'
#' ## save model parameters for documentation or dissimination
#' ## parameters can be load()-ed and used to predict the taper
#' ## or volume using one or several measured dbh
#' #spruce.taper.pars <- taper.model$par.lme
#' #save(spruce.taper.pars, file="spruce.taper.pars.rdata")


TapeR_FIT_LME.f <-
function(Id, x, y, knt_x, ord_x, knt_z, ord_z, IdKOVb = "pdSymm", control = list(), ...){
  ## control values, taken from nlme:::lme.formula()
  controlvals <- lmeControl()
  if (!missing(control)) {
    controlvals[names(control)] <- control
  }
	## require(nlme)
  BS_x <- BSplines(knots = knt_x, ord = ord_x, der = 0, x = x)
  BS_z <- BSplines(knots = knt_z, ord = ord_z, der = 0, x = x)

  # last column is omitted in order to ensure Ey(x=1)=0
  X <- BS_x[1:nrow(BS_x),1:ncol(BS_x)-1,drop=F]
  Z <- BS_z[1:nrow(BS_z),1:ncol(BS_z)-1,drop=F]

#   ****************************************************************************
#		fitting LME Modell E[y|X,Z] = X*b_fix + Z*b_rnd + eps
#   ****************************************************************************

  if(IdKOVb == "pdDiag"){ # VAR[b_rnd] from empiric KOV[b_rnd]
  	fit.lme <- lme(y ~ X-1,	random = list(Id = pdDiag(~Z-1)),
  	               control = controlvals)
  } else {
  	fit.lme <- lme(y ~ X-1,	random = list(Id = pdSymm(~Z-1)),
  	               control = controlvals, method = "ML" )
  }
#   ****************************************************************************

	b_fix <- as.numeric(fixef(fit.lme))
	KOVb_fix <- vcov(fit.lme)

	if(IdKOVb == "pdDiag"){ # empirische KOVb_rnd
		VARb_rnd <- matrix(getVarCov(fit.lme), ncol=ncol(getVarCov(fit.lme)))
		KORb_rnd <- cor(ranef(fit.lme))
		KOVb_rnd <- sqrt(VARb_rnd)%*%KORb_rnd%*%sqrt(VARb_rnd)
	} else {
		KOVb_rnd <- matrix(getVarCov(fit.lme), ncol=ncol(getVarCov(fit.lme)))
		KORb_rnd <- cov2cor(KOVb_rnd)
	}

	theta <- attr(fit.lme$apVar, "Pars")
	KOV_theta <- matrix(fit.lme$apVar, ncol=length(theta), byrow=T)

	sig2_eps <- as.numeric(exp(theta["lSigma"])^2)
  
	## calculate ratio between diameter in 1% and 5% of tree height
	## as an indicator about population average stump from
	f <- prod((XZ_BSPLINE.f(c(0.01, 0.05), knt_x, ord_x) %*% b_fix)^c(1,-1))

	par.lme <- list(knt_x = knt_x, pad_knt_x = TransKnots(knots=knt_x, ord=ord_x),
	                knt_z = knt_z, pad_knt_z = TransKnots(knots=knt_z, ord=ord_z),
	                ord_x = ord_x, ord_z = ord_z,
	                b_fix = b_fix, KOVb_fix = KOVb_fix,
	                sig2_eps = sig2_eps, dfRes = anova(fit.lme)$denDF,
	                KOVb_rnd = KOVb_rnd, theta = theta, KOV_theta = KOV_theta,
	                q001=f)

	return(list(fit.lme = fit.lme, par.lme = par.lme))
}
