#' @title qD.rout.f
#' @description Internal function not usually called by users
#' @param qD  vector of quantiles, finally passed to \code{pnorm}
#' @param alpha quantile for which root is sought
#' @param Hx Numeric vector of stem heights (m) along which to return the 
#' expected diameter
#' @param Hm Numeric vector of stem heights (m) along which diameter measurements 
#' were taken for calibration. Can be of length 1. Must be of same length as \code{Dm}
#' @param Dm Numeric vector of diameter measurements (cm) taken for calibration.
#'  Can be of length 1. Must be of same length as \code{Hm}
#' @param mHt Scalar. Tree height (m)
#' @param sHt Scalar. Standard deviation of stem height. Can be 0 if height was 
#' measured without error
#' @param par.lme List of taper model parameters obtained by 
#' \code{\link{TapeR_FIT_LME.f}}.
#' @param Rfn list with function name to provide estimated or assumed residual 
#' variances for the given measurements, optionally parameters for such functions
#' @param nGL Numeric scalar. Number of support points for numerical integration.
#' @param ... not currently used
#'
#' @return \code{qD} for given \code{alpha} with respect to \code{Int_CdN_DHx_dHt.f}
#' @author Edgar Kublin

qD.rout.f <-
function(qD, alpha = 0.975, Hx, Hm, Dm, mHt, sHt, par.lme, Rfn, nGL = 51, ...){

		Int_CdN_DHx_dHt.f(qD, Hx, Hm, Dm, mHt, sHt, par.lme, Rfn, nGL) - alpha
	}
