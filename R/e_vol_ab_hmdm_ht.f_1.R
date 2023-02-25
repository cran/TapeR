#' @title Estimate volume for stem and sections
#' @description Estimate volume for a complete stem from bottom to tip or 
#' for a section defined by lower and upper diameter or height. Variances for
#' estimated volumes are calculated.
#' @param Hm Numeric vector of stem heights (m) along which diameter 
#' measurements were taken for calibration. Can be of length 1. Must be of same 
#' length as \code{Dm}.
#' @param Dm Numeric vector of diameter measurements (cm) taken for calibration.
#' Can be of length 1. Must be of same length as \code{Hm}.
#' @param mHt Scalar. Tree height (m).
#' @param sHt Scalar. Standard deviation of stem height. Can be 0 if height was 
#' measured without error.
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
#' measured points (e.g. by \code{Rfn=list(fn="zero")}, cf. \code{\link{resVar}}).
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
#' @export
#'
#' @examples
#' # example data
#' data(DxHx.df)
#' # taper curve parameters based on all measured trees
#' data(SK.par.lme)
#' 
#' #select data of first tree
#' Idi <- (DxHx.df[,"Id"] == unique(DxHx.df$Id)[1])
#' (tree1 <- DxHx.df[Idi,])
#' 
#' ## Calculate the timber volume for the whole stem
#' VOL <- E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], 
#'                           Dm=tree1$Dx[3], 
#'                           mHt = tree1$Ht[1],
#'                           sHt = 0, # no height variance assumed
#'                           par.lme = SK.par.lme)
#' VOL$E_VOL #' expected value
#' VOL$VAR_VOL #' corresponding variance
#' VOL$Rfn
#' 
#' ## Calculate the timber volume for the whole stem, using Rfn="zero"
#' VOL <- E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], 
#'                           Dm=tree1$Dx[3], 
#'                           mHt = tree1$Ht[1],
#'                           sHt = 0, # no height variance assumed
#'                           par.lme = SK.par.lme,
#'                           Rfn = list(fn="zero"))
#' VOL$E_VOL #' expected value
#' VOL$VAR_VOL #' corresponding variance
#' VOL$Rfn
#' 
#' ## Calculate the timber volume for the whole stem
#' VOL <- E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], 
#'                           Dm=tree1$Dx[3], 
#'                           mHt = tree1$Ht[1],
#'                           sHt = 1, # no height variance assumed
#'                           par.lme = SK.par.lme)
#' VOL$E_VOL #' expected value
#' VOL$VAR_VOL #' corresponding variance
#' VOL$Rfn
#' 
#' ## Calculate the timber volume for the whole stem, using Rfn="zero"
#' VOL <- E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], 
#'                           Dm=tree1$Dx[3], 
#'                           mHt = tree1$Ht[1],
#'                           sHt = 1, #  height variance assumed
#'                           par.lme = SK.par.lme,
#'                           Rfn = list(fn="zero"))
#' VOL$E_VOL #' expected value
#' VOL$VAR_VOL #' corresponding variance
#' VOL$Rfn
#' 
#' ## Calculate the timber volume for a selected section given a height (0.3 - 5 m)
#' VOL <- E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], 
#'                           Dm=tree1$Dx[3], 
#'                           mHt = tree1$Ht[1], 
#'                           sHt = 1, 
#'                           par.lme = SK.par.lme, 
#'                           A=0.3, 
#'                           B=5, 
#'                           iDH = "H")
#' VOL$E_VOL #' expected value
#' VOL$VAR_VOL #' corresponding variance
#' VOL$Rfn
#' 
#' ## Calculate the timber volume for a selected section given a height (0.3 - 5 m)
#' VOL <- E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], 
#'                           Dm=tree1$Dx[3], 
#'                           mHt = tree1$Ht[1], 
#'                           sHt = 1, 
#'                           par.lme = SK.par.lme, 
#'                           A=0.3, 
#'                           B=5, 
#'                           iDH = "H", 
#'                           Rfn=list(fn="zero"))
#' VOL$E_VOL #' expected value
#' VOL$VAR_VOL #' corresponding variance
#' VOL$Rfn
#' 
#' ## Calculate the timber volume for a selected section given a diameter
#' ## threshold (30cm - 15cm) (negative value if A<B)
#' VOL <- E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3],
#'                           Dm=tree1$Dx[3], 
#'                           mHt = tree1$Ht[1], 
#'                           sHt = 1, 
#'                           par.lme = SK.par.lme, 
#'                           A=30, 
#'                           B=15, 
#'                           iDH = "D")
#' VOL$E_VOL #' expected value
#' VOL$VAR_VOL #' corresponding variance
#' 
#' ## Not run: 
#' \donttest{
#' ## The variance estimate resulting from the tree height uncertainty using
#' ## a Normal approximation takes much longer...
#' ptm <- proc.time()
#' E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1], 
#'                    sHt = 1, par.lme = SK.par.lme, IA=FALSE)
#' proc.time() - ptm
#' 
#' 
#' ##... than the calculation using a 2-point distribution...
#' ptm <- proc.time()
#' E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1],
#'                    sHt = 1, par.lme = SK.par.lme, IA=TRUE)
#' proc.time() - ptm
#' 
#' ##...fastest if no height variance is assumed
#' ptm <- proc.time()
#' E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1],
#'                    sHt = 0, par.lme = SK.par.lme, IA=FALSE)
#' proc.time() - ptm
#' 
#' ## Also the number of supportive points for the numerical integration
#' ## influences the calculation time
#' ptm <- proc.time()
#' E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1],
#'                    sHt = 0, par.lme = SK.par.lme, IA=FALSE, nGL=10)
#' proc.time() - ptm
#' ##' End(Not run)
#' }

E_VOL_AB_HmDm_HT.f <-
function(Hm, Dm, mHt, sHt = 0, A = NULL, B = NULL, iDH = "D", par.lme, Rfn = list(fn="sig2"), IA = F, nGL = 51, ...){
#   ************************************************************************************************

#		Hm; Dm; mHt = mw_HtT; sHt = sd_HtT; a = NULL; b = 7					; iDH = "DH"; par.lme = SK.par.lme; IA = F; nGL = 51
#		Hm; Dm; mHt = mw_HtT; sHt = sd_HtT; a = NULL; b = Int_E_VOL_dHt$Hb	; iDH = "H"; par.lme = SK.par.lme; IA = F; nGL = 51

#   a - unterer Grenzdurchmesser/ -hoehe (iDH = "H")
#		b - oberer Grenzdurchmesser/  -hoehe

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

		if(sHt > 0){#	Hoehentarifvarianz - Int{VOLab|(Hm,Dm),Ht]dHt}

			Ht = max(Hm,mHt)

		#   ************************************************************************
			Int_VOLab = Int_E_VOL_AB_HmDm_HT_dHt.f(Hm, Dm, A, B, iDH, mw_HtT=mHt, 
			                                       sd_HtT=sHt, par.lme, Rfn=Rfn, IA, nGL)
		#   ************************************************************************

			E_VOLab = Int_VOLab$E_VOL
			VAR_VOLab = Int_VOLab$VAR_VOL

		} else { # RotationsIntegral ueber die kalibrierte Schaftkurve E[D(Hx)|(Hm,Dm),Ht]


		#   ************************************************************************
			SK_VOLab = SK_VOLab_EBLUP_LME.f(xm, ym, a, b, Ht, par.lme, Rfn)
		#   ************************************************************************

			E_VOLab = SK_VOLab$VOL
			VAR_VOLab = SK_VOLab$VAR_VOL

		}

		Ht = max(Hm,mHt)

		Ha = a*Ht
		Hb = b*Ht

		Da = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = a, par.lme, Rfn)$yp
		Db = SK_EBLUP_LME.f(xm = Hm/Ht, ym = Dm, xp = b, par.lme, Rfn)$yp

		return(list(E_VOL = E_VOLab, VAR_VOL = VAR_VOLab, Hm = Hm, Dm = Dm, Ht = Ht, 
		            Da = Da, Db = Db, Ha = a*Ht, Hb = b*Ht, Rfn=Rfn))
	}
