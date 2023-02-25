#' @title Int_E_VOL_AB_HmDm_HT_dHt.f
#' @description Internal function not usually called by users
#' @param Hm Numeric vector of stem heights (m) along which diameter measurements 
#' were taken for calibration. Can be of length 1. Must be of same length as \code{Dm}
#' @param Dm Numeric vector of diameter measurements (cm) taken for calibration.
#'  Can be of length 1. Must be of same length as \code{Hm}
#' @param A Numeric scalar defining the lower threshold of a stem section for volume
#' estimation. Depends on \code{iDH}. If \code{iDH} = "D", a diameter
#' (cm), if \code{iDH} = "H", a height (m). If NULL, section starts at lowest point.
#' @param B Numeric scalar defining the upper threshold of a stem section for volume
#' estimation. Depends on \code{iDH}. If \code{iDH} = "D", a diameter
#' (cm), if \code{iDH} = "H", a height (m). If NULL, section ends at tip.
#' @param iDH Character scalar. Either "D" or "H". Type of threshold for
#' section volume estimation. See \code{A} or \code{B}.
#' @param mw_HtT Scalar. Tree height (m)
#' @param sd_HtT Scalar. Standard deviation of stem height. Can be 0 if height was 
#' measured without error
#' @param par.lme List of taper model parameters obtained by \code{\link{TapeR_FIT_LME.f}} 
#' @param Rfn list with function name to provide estimated or assumed residual 
#' variances for the given measurements, optionally parameters for such functions
#' @param IA Logic scalar. If TRUE, variance calculation of height
#' estimate based on 2-point distribution. If FALSE, variance calculation of height
#' estimate based on Normal approximation.
#' @param nGL Numeric scalar. Number of support points for numerical integration.
#' @param ... not currently used
#' @details integrating the taper curve considering uncertainty of height 
#' measurement
#' @return list with expected volume, variance of volume and squared expected value
#' incorporating the uncertainty of height measurement
#' @author Edgar Kublin
#' @import pracma

Int_E_VOL_AB_HmDm_HT_dHt.f <-
function(Hm, Dm, A=NULL, B=NULL, iDH="D", mw_HtT, sd_HtT, par.lme, Rfn=list(fn="sig2"), IA=FALSE, nGL=51, ...){

		if(IA){	# Two Point Approximation Lappi(2006)

			ncc = 2
			cc = list(x = c(mw_HtT - sd_HtT,mw_HtT + sd_HtT),w = c(1,1))

			E_VOLab = E2_VOLab = VAR_VOLab = dN_Ht = rep(0,ncc);
			Int_E_VOLab 	= Int_E2_VOLab = Int_VAR_VOLab = 0

			for (i in 1:ncc){

			#   ------------------------------------------------------------------------------------
				VOL = E_VOL_AB_HmDm_Ht.f(Hm, Dm, mHt=cc$x[i], A, B, iDH, par.lme, Rfn)
			#   ------------------------------------------------------------------------------------

				E_VOLab[i] 		= as.numeric(VOL$E_VOL)
				E2_VOLab[i] 	= as.numeric(VOL$E_VOL)^2
				VAR_VOLab[i] 	= as.numeric(VOL$VAR_VOL)

				dN_Ht[i]  		= 0.5                           			#   ZweiPunktVerteilung

				Int_E_VOLab		= Int_E_VOLab+cc$w[i]*dN_Ht[i]*E_VOLab[i]
				Int_E2_VOLab  	= Int_E2_VOLab+cc$w[i]*dN_Ht[i]*E2_VOLab[i]
				Int_VAR_VOLab 	= Int_VAR_VOLab+cc$w[i]*dN_Ht[i]*VAR_VOLab[i]
			}

		}else{ # Numerische Integration (Gauss - Legendre) ueber die Hoehenverteilung

			ncc = nGL

			cca	= mw_HtT - 5*sd_HtT;
			ccb	= mw_HtT + 5*sd_HtT; 		#	pnorm(q = b, mean = mw_HtT, sd = sd_HtT, lower.tail = T)

			cc 	= gaussLegendre(ncc,cca,ccb)	#;	cc; #	ncc = length(cc$x); # gaussLegendre(3,-3,3)

			E_VOLab = E2_VOLab = VAR_VOLab = dN_Ht = rep(0,ncc);

			Int_E_VOLab = Int_E2_VOLab = Int_VAR_VOLab = 0

			for (i in 1:ncc){

		#       Ht[i] = cc$x[i]

			#   ------------------------------------------------------------------------------------
				VOL = E_VOL_AB_HmDm_Ht.f(Hm, Dm, mHt=cc$x[i], A, B, iDH, par.lme, Rfn)
			#   ------------------------------------------------------------------------------------

				E_VOLab[i] 		= as.numeric(VOL$E_VOL)
				E2_VOLab[i] 	= as.numeric(VOL$E_VOL)^2
				VAR_VOLab[i] 	= as.numeric(VOL$VAR_VOL)

				dN_Ht[i]  		= dN.f(x = cc$x[i], mw = mw_HtT, sd = sd_HtT)

				Int_E_VOLab	= Int_E_VOLab+cc$w[i]*dN_Ht[i]*E_VOLab[i]
				Int_E2_VOLab  = Int_E2_VOLab+cc$w[i]*dN_Ht[i]*E2_VOLab[i]
				Int_VAR_VOLab = Int_VAR_VOLab+cc$w[i]*dN_Ht[i]*VAR_VOLab[i]

			}
		}

		E_VOL   = Int_E_VOLab
		VAR_VOL = Int_VAR_VOLab + Int_E2_VOLab - Int_E_VOLab^2

		return(list(E_VOL = E_VOL, VAR_VOL = VAR_VOL, E2_VOL = Int_E2_VOLab))

	}
