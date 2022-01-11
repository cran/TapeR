#' @title Estimate volume for stem and sections
#' @description Internal function not usually called by users
#' @param Hm Numeric vector of stem heights (m) along which diameter 
#' measurements were taken for calibration. Can be of length 1. Must be of same 
#' length as \code{Dm}.
#' @param Dm Numeric vector of diameter measurements (cm) taken for calibration.
#' Can be of length 1. Must be of same length as \code{Hm}.
#' @param mHt Scalar. Tree height (m).
#' @param A Numeric scalar defining the lower threshold of a stem section for 
#' volume estimation. Depends on \code{iDH}. If \code{iDH} = "D", a diameter 
#' (cm), if \code{iDH} = "H", a height (m). If NULL, section starts at lowest 
#' point.
#' @param B Numeric scalar defining the upper threshold of a stem section for 
#' volume estimation. Depends on \code{iDH}. If \code{iDH} = "D", a diameter
#' (cm), if \code{iDH} = "H", a height (m). If NULL, section ends at tip.
#' @param iDH Character scalar. Either "D" or "H". Type of threshold for section
#' volume estimation. See \code{A} or \code{B}.
#' @param par.lme List of taper model parameters obtained by 
#' \code{\link{TapeR_FIT_LME.f}}.
#' @param Rfn list with function name to provide estimated or assumed residual 
#' variances for the given measurements, optionally parameters for such functions
#' @param IA Logic scalar. If TRUE, variance calculation of height estimate 
#' based on 2-point distribution. If FALSE, variance calculation of height
#' estimate based on Normal approximation.
#' @param nGL Numeric scalar. Number of support points for numerical
#' integration.
#' @param ... not currently used
#' @details calculates the volume for a complete stem or sections defined by
#' \code{A} and \code{B}, which might be defined as diameter or height. The
#' parameter \code{Rfn} can be used to force the taper curve through the 
#' measured points (if \code{Rfn=list(fn="zero")}).
#' This function is used inside the two-point-approximation and numerical 
#' integration of the uncertainty. Evaluates the estimated taper curve 
#' repeatedly for all potential heights according to height uncertainty.
#' @return a list holding nine elements:
#' \itemize{
#'  \item{E_VOL: }{Estimated volume (m^3).}
#'  \item{VAR_VOL: }{Variance of the volume estimate.}
#'  \item{Hm: }{Height of diameter measurement (m).}
#'  \item{Dm: }{Diameter measurement (cm).}
#'  \item{Ht: }{Tree height (m).}
#'  \item{Da: }{Diameter at lower section threshold (cm).}
#'  \item{Db: }{Diameter at upper section threshold (cm).}
#'  \item{Ha: }{Height at lower section threshold (m).}
#'  \item{Hb: }{Height at upper section threshold (m).}
#'  \item{Rfn: }{Function applied for estimated or assumed residual variance.}
#' }
#' @author Edgar Kublin
#' @references Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem 
#' taper and volume prediction method based on mixed-effects B-spline 
#' regression, Eur J For Res, 132:983-997.
#' @seealso \code{\link{TapeR_FIT_LME.f}}

E_VOL_AB_HmDm_Ht.f <-
function(Hm, Dm, mHt, A = NULL, B = NULL, iDH = "D", par.lme, Rfn=list(fn="sig2"), ...){

#   A - unterer Grenzdurchmesser/ -hoehe
#		B - oberer Grenzdurchmesser / -hoehe

		Ht = max(Hm,mHt)

		if(min(Dm)>0){Ht = max(c(Hm,Ht))}else{Ht = max(Hm)}

		xm = Hm/Ht
		ym = Dm

		if(is.null(A)){
			a=0
		}else{
			if(iDH %in% c("d","D")){
				a = xy0_SK_EBLUP_LME.f(xm, ym, y0 = A, par.lme, Rfn)
			}else{
				a = min(1,A/Ht)
			}
		}

		if(is.null(B)){
			b=1
		}else{
			if(iDH %in% c("d","D")){
				b = xy0_SK_EBLUP_LME.f(xm, ym, y0 = B, par.lme, Rfn)
			}else{
				b = min(1,B/Ht)
			}
		}

#       Abschnittsvolumen zu den Kalibrierungsdaten (Hm,Dm) und Schafthoehe Ht 

#       ----------------------------------------------------------
		SK_VOLab = SK_VOLab_EBLUP_LME.f(xm, ym, a, b, Ht, par.lme, Rfn = Rfn)
#       ----------------------------------------------------------

		E_VOLab = SK_VOLab$VOL
		VAR_VOLab = SK_VOLab$VAR_VOL

		Ht = max(Hm, mHt)

		Ha = a*Ht
		Hb = b*Ht

		Da = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = a, par.lme, Rfn)$yp
		Db = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = b, par.lme, Rfn)$yp

		return(list(E_VOL = E_VOLab,VAR_VOL = VAR_VOLab, Hm = Hm, Dm = Dm, 
		            Ht = Ht, Da = Da, Db = Db, Ha = a*Ht, Hb = b*Ht, Rfn=Rfn))

	}
