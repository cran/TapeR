#' @title Int_CdN_DHx_dHt.f
#' @description Internal function not usually called by users
#' @param qD vector of quantiles, finally passed to \code{pnorm}
#' @param Hx Numeric vector of stem heights (m) along which to return the 
#' expected diameter
#' @param Hm Numeric vector of stem heights (m) along which diameter measurements 
#' were taken for calibration. Can be of length 1. Must be of same length as \code{Dm}
#' @param Dm Numeric vector of diameter measurements (cm) taken for calibration.
#'  Can be of length 1. Must be of same length as \code{Hm}
#' @param mHt Scalar. Tree height (m)
#' @param sHt Scalar. Standard deviation of stem height. Can be 0 if height was 
#' measured without error
#' @param par.lme List of taper model parameters obtained by \code{\link{TapeR_FIT_LME.f}}
#' @param Rfn list with function name to provide estimated or assumed residual 
#' variances for the given measurements, optionally parameters for such functions
#' @param nGL Numeric scalar. Number of support points for numerical integration
#' @param ... not currently used
#'
#' @return Int_CdN_dN
#' @author Edgar Kublin
#' @import pracma

Int_CdN_DHx_dHt.f <-
function(qD, Hx, Hm, Dm, mHt, sHt, par.lme, Rfn, nGL = 51, ...){
#   ------------------------------------------------------------------------------------------------
	# Hx=Hx[i]
		ncc = nGL

		cca	= mHt - 5*sHt;
		ccb	= mHt + 5*sHt; 		#	pnorm(q = b, mean = mw_HtT, sd = sd_HtT, lower.tail = T)

		cc 	= gaussLegendre(ncc,cca,ccb)	#;	cc; #	ncc = length(cc$x); # gaussLegendre(3,-3,3)

		Mw_DHxHt = StD_DHxHt = dN_Ht = CdN_DHxHt = w_CdN_dN = Sum_w_CdN_dN =rep(0,ncc);

		Int_CdN_dN = 0

		for (i in 1:ncc){

	#       Ht[i] = cc$x[i]

			SK 	= E_DHx_HmDm_HT.f( Hx, Hm, Dm, mHt = cc$x[i], sHt = 0, par.lme, Rfn)

			Mw_DHxHt[i] 	= as.numeric(SK$DHx)
			StD_DHxHt[i] 	= sqrt(as.numeric(SK$MSE_Mean))

			dN_Ht[i]  		= dN.f(x = cc$x[i], mw = mHt, sd = sHt)
			CdN_DHxHt[i] 	= CdN_DHxHt.f(Ht = cc$x[i], Hx, qD, Hm, Dm, par.lme, Rfn)

			w_CdN_dN[i]		= cc$w[i]*dN_Ht[i]*CdN_DHxHt[i]

			Sum_w_CdN_dN[i] = Int_CdN_dN + w_CdN_dN[i]

			Int_CdN_dN = Sum_w_CdN_dN[i]
		}

		return(Int_CdN_dN)

	}
