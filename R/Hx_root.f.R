#' @title find root (height) given diameter, measurements and fitted model
#' @description Internal function not usually called by users
#' @param Hx Numeric vector of stem heights (m) along which to return the 
#' expected diameter
#' @param Dx expected diameter
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
#' @param ... not currently used
#' @details function is called by \code{\link[stats]{uniroot}} inside 
#' \code{\link{E_HDx_HmDm_HT.f}}
#' @return deviation between observed diameter \code{Dx} and diameter in height
#' \code{Hx}.
#' @author Edgar Kublin

Hx_root.f <- function(Hx,Dx,Hm,Dm,mHt,sHt,par.lme, Rfn, ...){

		return(E_DHx_HmDm_HT.f( Hx, Hm, Dm, mHt, sHt = 0, par.lme, Rfn)$DHx - Dx)

	}
