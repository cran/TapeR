#' @title find relative height of given diameter
#' @description Internal function not usually called by users
#' @param xm relative heights for which measurements are available
#' @param ym corresponding diameter measurements in height \code{xm}
#' @param y0 given diameter for which height is required
#' @param par.lme Fitted model object, return of \code{\link{TapeR_FIT_LME.f}}
#' @param Rfn list with function name to provide estimated or assumed residual 
#' variances for the given measurements, optionally parameters for such functions
#' @param ... not currently used
#' @details function used to transform given diameter in volume calculation into
#' height; c.f \code{\link{E_VOL_AB_HmDm_HT.f}}; with \code{Rfn} one can decide 
#' whether the measured diameters are forced to lie exactly on the taper curve
#' \code{Rfn$fn="sig2"} or not \code{Rfn$fn="zero"}. Other options are possible, 
#' see also \code{\link{SK_EBLUP_LME.f}} and \code{\link{resVar}}.
#' @return relative height of given diameter \code{y0}
#' @author Edgar Kublin
#' @importFrom stats uniroot

xy0_SK_EBLUP_LME.f <-
function(xm, ym, y0, par.lme, Rfn=list(fn="sig2"), ...){

		SK_LME = SK_EBLUP_LME.f(xm = xm, ym = ym, xp = c(0), 
		                        par.lme = par.lme, Rfn=Rfn)

		xmin <- uniroot(xy0_root.f, c(0,1), tol = 0.00001, y0 = y0, 
		                SK = SK_LME, par.lme = par.lme)

		x0 = xmin$root

		return(x0)

	}
