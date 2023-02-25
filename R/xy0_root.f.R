#' @title xy0_root.f
#' @description Internal function not usually called by users
#' @param x relative height
#' @param y0 diameter for which height is required
#' @param SK return of \code{\link{SK_EBLUP_LME.f}} containing estimated fixed
#'  and the tree specific random effects of the taper model
#' @param par.lme List of taper model parameters obtained by \code{\link{TapeR_FIT_LME.f}} 
#' @param ... not currently used
#' @details used in \code{xy0_SK_EBLUP_LME.f} to find the root of taper curve
#' (i.e. height \code{x}) at given diameter \code{y0}
#'
#' @return difference between actual diameter at height \code{x} and given
#' diameter \code{y0}
#' @author Edgar Kublin

xy0_root.f <-
function (x, y0, SK, par.lme, ...){

		b_fix = SK$b_fix
		b_rnd = SK$b_rnd

		X = XZ_BSPLINE.f(x, par.lme$knt_x, par.lme$ord_x)
		Z = XZ_BSPLINE.f(x, par.lme$knt_z, par.lme$ord_z)

		xy0_root = X%*%b_fix + Z%*%b_rnd - y0

		return(xy0_root)
	}
