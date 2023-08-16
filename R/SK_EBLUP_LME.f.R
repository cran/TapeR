#' @title Evaluate fitted taper curve
#' @description This is the actual function to estimate diameters according to
#' the fitted mixed B-splines model.
#' @param xm relative heights for which measurements are available
#' @param ym corresponding diameter measurements in height \code{xm}
#' @param xp relative heights for which predictions are required
#' @param par.lme Fitted model object, return of \code{\link{TapeR_FIT_LME.f}}
#' @param Rfn list with function name to provide estimated or assumed residual 
#' variances for the given measurements, optionally parameters for such functions
#' @param ... not currently used
#' @details This function is the actual working horse for prediction using the
#' fitted taper model. Based on the model \code{par.lme} and the measured
#' diameters \code{ym} and corresponding (relative) heights \code{xm} of a 
#' specific tree (there might be just one measurement), the random
#' effect parameters and subsequently diameters are estimated. Depending on the
#' parameter \code{Rfn}, the calibrated taper curve is forced through the
#' given diameter \code{ym} (\code{Rfn = list(fn="zero")}), or calibrated using
#' the complete residual variance-covariance information 
#' (\code{Rfn = list(fn="sig2")}, the default). 
#' Further assumptions are possible, see also \code{\link{resVar}} and
#' Kublin et al. (2013) p. 987 for more details.
#' @return a list holding nine elements:
#' \itemize{
#'   \item b_fix fixed effects parameter of taper model
#'   \item b_rnd random effects parameter given tree (posterior mean b_k)
#'   \item yp estimated diameter in height \code{xp}
#'   \item KOV_Mean variance-covariance matrix of expected value
#'   \item KOV_Pred variance-covariance matrix of prediction
#'   \item CI_Mean mean and limits of confidence interval
#'   \item MSE_Mean mean squared error of expected value
#'   \item MSE_Pred mean squared error of prediction
#'   \item CI_Pred mean and limits of prediction interval
#' }
#' @author Edgar Kublin
#' @seealso \code{\link{E_DHx_HmDm_HT.f}}, \code{\link{E_VOL_AB_HmDm_HT.f}}, 
#' \code{\link{resVar}}
#' @importFrom stats qt
#' @examples
#' data("SK.par.lme")
#' TapeR:::SK_EBLUP_LME.f(1.3/27, 30, 1.3/27, SK.par.lme)
#' ## using empirical best linear unbiased estimator: estimate != 30
#' TapeR:::SK_EBLUP_LME.f(1.3/27, 30, 1.3/27, SK.par.lme, Rfn=list(fn="sig2"))$yp
#' ## interpolate / force through given diameter: estimate  == 30
#' TapeR:::SK_EBLUP_LME.f(1.3/27, 30, 1.3/27, SK.par.lme, Rfn=list(fn="zero"))$yp
#' TapeR:::SK_EBLUP_LME.f(1.3/27, 30, c(1.3, 5)/27, SK.par.lme)
#' par.lme <- SK.par.lme
#' h <- 12 # tree height
#' xm <- c(1.3, 3) / h # relative measuring height
#' ym <- c(8, 7.5) # measured diameter
#' xp <- c(0.5, 1) / h # relative prediction height
#' TapeR:::SK_EBLUP_LME.f(xm, ym, xp, SK.par.lme)
#' @export

SK_EBLUP_LME.f <-
function(xm, ym, xp, par.lme, Rfn=list(fn="sig2"), ...){
#   ************************************************************************************************

	#	xm = xm_i; ym = ym_i; xp = xp_i; par.lme = SK_FIT_LME$par.lme

		EBLUP_b_k = MSE_Mean = MSE_Pred = CI_Mean = CI_Pred = NULL

	#   Design Matrizen X und Z zu den Kalibrierungsdaten :.........................................

		x_k 	= xm
		y_k   = ym

		X_k = XZ_BSPLINE.f(x_k, par.lme$knt_x, par.lme$ord_x)
		Z_k = XZ_BSPLINE.f(x_k, par.lme$knt_z, par.lme$ord_z)

	#   Feste Effekte - Mittlere Schaftkurve in der GesamtPopulation (PA) : E[y|x] = X*b_fix:.......

		b_fix 		= par.lme$b_fix
		KOVb_fix    = par.lme$KOVb_fix

	#   Kovarianz fuer die Random Effekte (fit.lme):................................................

		KOV_b 		= par.lme$KOVb_rnd
		Z_KOVb_Zt_k = Z_k%*%KOV_b%*%t(Z_k)                                			# fix(Z_KOVb_Zt)

	#   Residuen Varianz:...........................................................................

		sig2_eps = par.lme$sig2_eps
		dfRes    = par.lme$dfRes
    
	  rv <- resVar(x_k, fn = Rfn$fn, sig2 = sig2_eps, par = Rfn$par)
	  R_k = diag(rv, ncol(Z_KOVb_Zt_k))
	  
	#   Kovarianzmatrix der Beobachtungen (Sigma.i) :...............................................

		KOV_y_k		= Z_KOVb_Zt_k + R_k                       	#   SIGMA_k in V&C(1997)(6.2.3)
		KOVinv_y_k  = solve(KOV_y_k);                    		#   SIGMA_k^(-1)

	#   EBLUP - Posterior mean (b_k) aus Einhaengung (xm,ym) berechnen :.............................

	#   ***************************************************************
		EBLUP_b_k 	= KOV_b%*%t(Z_k)%*%KOVinv_y_k%*%(y_k - X_k%*%b_fix); #   V&C(1997) (6.2.49)
	#   ***************************************************************

	#	x_pre 	= unique(xp[order(xp)])
	#	x_pre 	= xp[order(xp)]
		x_pre 	= xp

		X_kh = XZ_BSPLINE.f(x_pre, par.lme$knt_x, par.lme$ord_x)
		Z_kh = XZ_BSPLINE.f(x_pre, par.lme$knt_z, par.lme$ord_z)

	#   ********************************************************************************************
		yp = EBLUP_y_kh  = X_kh%*%b_fix + Z_kh%*%EBLUP_b_k    		#	V&C(1997) (6.2.52)
	#   ********************************************************************************************

	#   ----------------------------------------------------------------------------------------
	#   		Vorhersageintervalle Schaftkurve lmeBLUP (SS) mit Einhaengung in (x_k,y_k)
	#   ----------------------------------------------------------------------------------------

	if(T){

	#	rm(v_k, V_k, C_k)

	#   Posterior Varianz V^*_k = VAR[b_k|y_k,beta,tetha(sig2_eps)] (V&C s.252 6.2.50) :............

		Vv_k = KOV_b - KOV_b%*%t(Z_k)%*%KOVinv_y_k%*%Z_k%*%KOV_b

	#   Vorhersage Varianz: KOV[(^b_k-b_k)|y_k,beta,tetha(sig2_eps)] - V&C (6.2.51):...................

		V_k = Vv_k + KOV_b%*%t(Z_k)%*%KOVinv_y_k%*%X_k%*%KOVb_fix%*%t(X_k)%*%KOVinv_y_k%*%Z_k%*%KOV_b

	#   KOV[(beta^-beta),(b^_k-b_k)] V&C (1997) (6.2.54):...........................................

		C_k = -KOVb_fix%*%t(X_k)%*%KOVinv_y_k%*%Z_k%*%KOV_b

	#   MSE(Mittelwert/Vorhersage) :................................................................

	#	rm(KOV_Mean,KOV_Pred,MSE_Mean,MSE_Pred,R_kh)

		R_kh = diag(sig2_eps,nrow(X_kh))					#   Anzahl Beobachtungen

	#   ********************************************************************************************
		KOV_Mean = X_kh%*%KOVb_fix%*%t(X_kh) + Z_kh%*%V_k%*%t(Z_kh) + X_kh%*%C_k%*%t(Z_kh) + Z_kh%*%t(C_k)%*%t(X_kh)
		KOV_Pred = KOV_Mean + R_kh      							#   V&C(1997) (6.2.55)
	#   ********************************************************************************************

	#	KOV_Mean = round(KOV_Mean, digits=4)

		MSE_Mean = round(diag(KOV_Mean),digits=4)
		MSE_Pred = round(diag(KOV_Pred),digits=4)

		c_alpha = qt(p=0.025, df=dfRes, ncp=0, lower.tail = F, log.p = FALSE)

		CI_Mean = cbind(EBLUP_y_kh - c_alpha*sqrt(MSE_Mean),EBLUP_y_kh,EBLUP_y_kh + c_alpha*sqrt(MSE_Mean))
		CI_Pred = cbind(EBLUP_y_kh - c_alpha*sqrt(MSE_Pred),EBLUP_y_kh,EBLUP_y_kh + c_alpha*sqrt(MSE_Pred))
	}

	#   ********************************************************************************************
		return(list(b_fix 	= b_fix, 		b_rnd 	 = EBLUP_b_k,
					yp 		= EBLUP_y_kh, 	KOV_Mean = KOV_Mean, KOV_Pred = KOV_Pred, CI_Mean = CI_Mean,
											MSE_Mean = MSE_Mean, MSE_Pred = MSE_Pred, CI_Pred = CI_Pred
					))
	#   ********************************************************************************************

	}
