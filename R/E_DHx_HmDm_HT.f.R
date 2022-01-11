#' @title Estimate diameter and approximate confidence and prediction intervals
#' @description Calibrates a taper curve based on at least one diameter 
#' measurement and returns the expected diameters and approximate variances
#' @param Hx Numeric vector of stem heights (m) along which to return the 
#' expected diameter.
#' @param Hm Numeric vector of stem heights (m) along which diameter 
#' measurements were taken for calibration. Can be of length 1. Must be of same 
#' length as \code{Dm}.
#' @param Dm Numeric vector of diameter measurements (cm) taken for calibration.
#' Can be of length 1. Must be of same length as \code{Hm}.
#' @param mHt Scalar. Tree height (m).
#' @param sHt Scalar. Standard deviation of stem height. Can be 0 if height was 
#' measured without error.
#' @param par.lme List of taper model parameters obtained by 
#' \code{\link{TapeR_FIT_LME.f}}.
#' @param Rfn list with function name to provide estimated or assumed residual 
#' variances for the given measurements, optionally parameters for such functions
#' @param ... not currently used
#' @details calibrates the tree specific taper curve and calculates approximate
#' confidence intervals, which can be useful for plotting. Uncertainty resulting
#' from tariff height estimates if tree height was not measured is incorporated.
#' Using \code{Rfn} the taper curve can be forced through the measured 
#' diameters, c.f. \code{\link{resVar}}.
#' @return a list holding six elements:
#' \itemize{
#'  \item{DHx: }{Numeric vector of diameters (cm) (expected value) along the 
#'  heights given by \code{Hx}.}
#'  \item{Hx: }{Numeric vector of heights (m) along which to return the expected
#'   diameter.}
#'  \item{MSE_Mean: }{Mean squared error for the expected value of the diameter.}
#'  \item{CI_Mean: }{Confidence interval. Matrix of the 95\% conf. int. for the 
#'  expected value of the diameter (cm). First column: lower limit, second 
#'  column: mean, third column: upper limit.}
#'  \item{MSE_Pred: }{Mean squared error for the prediction of the diameter.}
#'  \item{CI_Mean: }{Prediction interval. Matrix of the 95\% conf. int. for the 
#'  prediction of the diameter (cm). First column: lower limit, second column: 
#'  mean, third column: upper limit.}
#'  \item{Rfn: }{Function applied for estimated or assumed residual variance.}
#' }
#' @author Edgar Kublin
#' @references Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem 
#' taper and volume prediction method based on mixed-effects B-spline 
#' regression, Eur J For Res, 132:983-997.
#' @seealso \code{\link{TapeR_FIT_LME.f}}
#' @importFrom stats qt predict
#' @importFrom graphics matlines
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
#' ## Predict the taper curve based on the diameter measurement in 2 m
#' ## height and known height 
#' tc.tree1 <- E_DHx_HmDm_HT.f(Hx=1:tree1$Ht[1], 
#'                             Hm=tree1$Hx[3],
#'                             Dm=tree1$Dx[3], 
#'                             mHt = tree1$Ht[1], 
#'                             sHt = 0, 
#'                             par.lme = SK.par.lme)
#' #plot the predicted taper curve
#' plot(tc.tree1$Hx, tc.tree1$DHx, type="l", las=1)
#' #lower CI
#' lines(tc.tree1$Hx, tc.tree1$CI_Mean[,1], lty=2)
#' #upper CI
#' lines(tc.tree1$Hx, tc.tree1$CI_Mean[,3], lty=2)
#' #lower prediction interval
#' lines(tc.tree1$Hx, tc.tree1$CI_Pred[,1], lty=3)
#' #upper prediction interval
#' lines(tc.tree1$Hx, tc.tree1$CI_Pred[,3], lty=3)
#' #add measured diameter used for calibration
#' points(tree1$Hx[3], tree1$Dx[3], pch=3, col=2)
#' #add the observations
#' points(tree1$Hx, tree1$Dx)
#' 
#' ## feature of forcing taper curve through measured diameters
#' i <- c(3, 6)
#' tc.tree1 <- E_DHx_HmDm_HT.f(Hx=seq(0, tree1$Ht[1], 0.1), 
#'                             Hm=tree1$Hx[i],
#'                             Dm=tree1$Dx[i], 
#'                             mHt = tree1$Ht[1], 
#'                             sHt = 0, 
#'                             par.lme = SK.par.lme,
#'                             Rfn=list(fn="sig2"))
#' tc.tree2 <- E_DHx_HmDm_HT.f(Hx=seq(0, tree1$Ht[1], 0.1), 
#'                             Hm=tree1$Hx[i],
#'                             Dm=tree1$Dx[i], 
#'                             mHt = tree1$Ht[1], 
#'                             sHt = 0, 
#'                             par.lme = SK.par.lme,
#'                             Rfn=list(fn="zero"))
#' # plot the predicted taper curve
#' plot(tc.tree1$Hx, tc.tree1$DHx, type="l", las=1)
#' # added taper curve through measurement
#' points(x=tc.tree2$Hx, y=tc.tree2$DHx, type="l", lty=2)
#' # closer window
#' plot(tc.tree1$Hx, tc.tree1$DHx, type="l", las=1, xlim=c(0, 8), ylim=c(24, 30))
#' # added taper curve through measurement
#' points(x=tc.tree2$Hx, y=tc.tree2$DHx, type="l", lty=2)
#' # add measured diameter used for calibration
#' points(tree1$Hx[i], tree1$Dx[i], pch=3, col=2)
#' # add the observations
#' points(tree1$Hx, tree1$Dx)
#' 
#' ## apply yet another residual variance function
#' i <- c(1, 2, 3) # calibrating with 0.5, 1m and 2m, assuming no error in 0.5m
#' zrv <- tree1$Hx[1] / tree1$Ht[1] # assumed zero resiudal variance
#' # assumed residual variance per measurement
#' TapeR:::resVar(relH = tree1$Hx[i] / tree1$Ht[1], fn = "dlnorm", 
#'                sig2 = SK.par.lme$sig2_eps, par = list(zrv=zrv))
#' tc.tree3 <- E_DHx_HmDm_HT.f(Hx=seq(0, tree1$Ht[1], 0.1), 
#'                             Hm=tree1$Hx[i],
#'                             Dm=tree1$Dx[i], 
#'                             mHt = tree1$Ht[1], 
#'                             sHt = 0, 
#'                             par.lme = SK.par.lme,
#'                             Rfn=list(fn="dlnorm", par=list(zrv=zrv)))
#' plot(tc.tree1$Hx, tc.tree1$DHx, type="l", las=1, xlim=c(0, 4))
#' points(x=tc.tree3$Hx, y=tc.tree3$DHx, type="l", lty=2)
#' points(tree1$Hx[i], tree1$Dx[i], pch=3, col=2)
#' points(tree1$Hx, tree1$Dx)

E_DHx_HmDm_HT.f <-
function( Hx, Hm, Dm, mHt, sHt = 0, par.lme, Rfn=list(fn="sig2"), ...){

		if(sHt > 0){

		#	SK.par.lme = par.lme

			mw_HtT 	= mHt
			sd_HtT 	= sHt

			nx 		= c(30,20)

			Ht_u	= mw_HtT - 2*sd_HtT
			Ht_m    = mw_HtT
			Ht_o	= mw_HtT + 2*sd_HtT

		#   ------------------------------------------------------------------------
		#   										CI_SK_u
		#   ------------------------------------------------------------------------

		#   SK_u :..................................................................

			hx      = seq(0,1,length.out = sum(nx))

			SK_u 	= SK_EBLUP_LME.f(xm = Hm/Ht_u, ym = Dm, xp = hx, par.lme, Rfn) # Kalibrierung/LME

			u.ssp   = smooth.spline(x = hx*Ht_u, y = SK_u$yp, df=10)
					  u.ssp$fit$coef[length(u.ssp$fit$coef)] = 0    			#   u.ssp(Ht_u)=0

		#   CI_u(SK_m) :............................................................

			SK_m    = SK_EBLUP_LME.f(xm = Hm/Ht_m, ym = Dm, xp = hx, par.lme, Rfn)

			m.ssp   = smooth.spline(x = hx*Ht_m,  y = SK_m$CI_Pred[,1] , df=10)
					  m.ssp$fit$coef[length(u.ssp$fit$coef)] = 0    			#   m.ssp(Ht_m)=0

		#   CI_SK_u (Approximation):................................................

		#	nx = c(30,20)

			hx          = c(seq(0,0.3*Ht_u,length.out = nx[1]),seq(0.6*Ht_u,Ht_u,length.out = nx[2]))

			qu        	= predict(u.ssp,x = hx, deriv = 0)$y;   	qu 		= apply(cbind(0,qu),1,max)
			qm        	= predict(m.ssp,x = hx, deriv = 0)$y;		qm 		= apply(cbind(0,qm),1,max)
			qmin        = apply(cbind(qu,qm),1,min);				qmin 	= apply(cbind(qmin,0),1,max)

		  qD_u.ssp	= smooth.spline(x = hx,  y = qmin, df=10); qD_u.ssp$fit$coef[length(qD_u.ssp$fit$coef)] = 0

		#   ----------------------------------------------------------------------------------------
		#   										CI_SK_o
		#   ----------------------------------------------------------------------------------------

			hx      	= seq(0,1,length.out = sum(nx))

			SK_o		= SK_EBLUP_LME.f(xm = Hm/Ht_o, ym = Dm, xp = hx, par.lme, Rfn)

			o.ssp       = smooth.spline(x = hx*Ht_o,  y = SK_o$yp, df=10); o.ssp$fit$coef[length(u.ssp$fit$coef)] = 0    #   u.ssp(Ht_u)=0

		#   ----------------------------------------------------------------------------------------

			SK_m       	= SK_EBLUP_LME.f(xm = Hm/Ht_m, ym = Dm, xp = hx, par.lme, Rfn)
			m.ssp       = smooth.spline(x = hx*Ht_m,  y = SK_m$CI_Pred[,3] , df=10)

		#   CI_SK_o (Approximation):................................................................

		#	nx=c(30,20)

			hx          = c(seq(0,0.3*Ht_o,length.out = nx[1]),seq(min(Ht_u,0.6*Ht_o),Ht_o,length.out=nx[2]))

			qm        	= predict(m.ssp,x = hx, deriv = 0)$y;	qm = apply(cbind(0,qm),1,max)
			qo        	= predict(o.ssp,x = hx, deriv = 0)$y;	qo = apply(cbind(0,qo),1,max)

			qmax        = apply(cbind(qm,qo),1,max)

		    qD_o.ssp	= smooth.spline(x = hx,  y = qmax, df=10); qD_o.ssp$fit$coef[length(qD_o.ssp$fit$coef)] = 0

		#   ----------------------------------------------------------------------------------------


			HHx = Hx[order(Hx)]

			qD_u = qD_o = list(x=0,y=0)

			hx = HHx/Ht_m; hx = apply(cbind(hx,1),1,min)

			SK_m   	 = SK_EBLUP_LME.f(xm = Hm/Ht_m, ym = Dm, xp = hx, par.lme, Rfn)

			DHx 	 = SK_m$yp
			MSE_Mean = SK_m$MSE_Mean
			MSE_Pred = SK_m$MSE_Pred

			HHx = Hx[order(Hx)]; HHx = apply(cbind(Ht_u,HHx),1,min)

			qD_u[["x"]] = HHx
			qD_u[["y"]] = predict(qD_u.ssp,x = HHx)$y

			HHx = Hx[order(Hx)]; HHx = apply(cbind(Ht_o,HHx),1,min)

			qD_o[["x"]] = HHx
			qD_o[["y"]] = predict(qD_o.ssp,x = HHx)$y

			CI_Mean = cbind(qD_u[["y"]],DHx,qD_o[["y"]])

			sig2_eps  = par.lme$sig2_eps
			dfRes    	= par.lme$dfRes
			c_alpha 	= qt(p=0.025, df=dfRes, ncp=0, lower.tail = F, log.p = FALSE)

			CI_Pred 	= CI_Mean

		#	qD_u = CI_Mean[,2] - c_alpha*sqrt(MSE_Mean)
		#   sqrt(MSE_Mean) = (CI_Mean[,2]-qD_u)/c_alpha

			dci      	= sqrt(((DHx - qD_u[["y"]])/c_alpha)^2 + sig2_eps)
			CI_Pred[,1] = DHx - c_alpha*dci
			CI_Pred[,1]	= apply(cbind(CI_Pred[,1],0),1,max)

		#	ii = (Ht_u<=Hx); CI_Pred[ii,1] = 0

			dci        	= sqrt(((qD_o[["y"]] - DHx)/c_alpha)^2 + sig2_eps)
			CI_Pred[,3] = DHx + c_alpha*dci
			CI_Pred[,3]	= apply(cbind(CI_Pred[,3],0),1,max)

		#	ii = (Ht_o<=Hx); CI_Pred[ii,3]=0

			if(F){
				plot(Hx,CI_Mean[,3],type="n")
				matlines(Hx,CI_Mean, col = "blue", lwd=2, lty=1)
				matlines(Hx,CI_Pred,col = "red", lwd=2, lty=1)

				cbind(CI_Pred[,1],CI_Mean[,1])
				cbind(CI_Mean[,3],CI_Pred[,3])
			}

		}else{ #  mHt=DxHx.df[Idi,"Ht"][1]gemessen

			xm = Hm/mHt
			ym = Dm

			hx = Hx[order(Hx)]/mHt
			# hx = apply(cbind(1,hx),1,min)
			# hx = hx[hx <= 1]
			hx = ifelse(hx > 1, 1, hx)

		#   ----------------------------------------------------------------------------------------
			SK_m = SK_EBLUP_LME.f(xm = Hm/mHt, ym = Dm, xp = hx, par.lme, Rfn) # Kalibrierung/LME
		#   ----------------------------------------------------------------------------------------

			DHx				= 	SK_m$yp

			MSE_Mean		= 	SK_m$MSE_Mean
			MSE_Pred		= 	SK_m$MSE_Pred

			dfRes    = par.lme$dfRes
			c_alpha = qt(p=0.025, df=dfRes, ncp=0, lower.tail = F, log.p = FALSE)

			CI_Mean = cbind(DHx - c_alpha*sqrt(MSE_Mean),DHx,DHx + c_alpha*sqrt(MSE_Mean))
			CI_Pred = cbind(DHx - c_alpha*sqrt(MSE_Pred),DHx,DHx + c_alpha*sqrt(MSE_Pred))

		}


		return(list(DHx = DHx, Hx = Hx[order(Hx)],
		            MSE_Mean = MSE_Mean, CI_Mean = CI_Mean,
		            MSE_Pred = MSE_Pred, CI_Pred = CI_Pred,
		            Rfn=Rfn))

  }
