#' @title Estimate diameter and exact confidence and prediction intervals
#' @description Calibrates a taper curve based on at least one diameter
#' measurement and returns the expected diameters and exact variances
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
#' @details calibrates the tree specific taper curve and calculates 'exact'
#' confidence intervals, which can be useful for plotting.
#' Attention: this function is somewhat time-consuming.
#' @return a matrix with six columns:
#' \itemize{
#'  \item{Hx: }{Numeric vector of heights (m) along which to return the expected
#'  diameter.}
#'  \item{q_DHx_u: }{Lower confidence interval (cm). (95\% CI except for estimates
#'  close to the stem tip.)}
#'  \item{DHx: }{Diameter estimate (cm).}
#'  \item{q_DHx_o: }{Upper CI (cm).}
#'  \item{cP_DHx_u: }{Probability of observations <\code{q_DHx_u}.}
#'  \item{cP_DHx_o: }{Probability of observations <\code{q_DHx_o}.}
#' }
#' @author Edgar Kublin
#' @references Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem
#' taper and volume prediction method based on mixed-effects B-spline
#' regression, Eur J For Res, 132:983-997.
#' @seealso \code{\link{TapeR_FIT_LME.f}}
#' @importFrom stats uniroot qt
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
#' ## Calculate "exact" CIs. Careful: This takes a while!
#' #library(pracma)# for numerical integration with gaussLegendre()
#' \donttest{
#' tc.tree1.exact <- E_DHx_HmDm_HT_CIdHt.f(Hx=1:tree1$Ht[1],
#'                                         Hm=tree1$Hx[3],
#'                                         Dm=tree1$Dx[3],
#'                                         mHt=tree1$Ht[1],
#'                                         sHt=1,
#'                                         par.lme=SK.par.lme)
#' #add exact confidence intervals to approximate intervals above - fits
#' #quite well
#' lines(tc.tree1.exact[,1], tc.tree1.exact[,2], lty=2,col=2)
#' lines(tc.tree1.exact[,1], tc.tree1.exact[,4], lty=2,col=2)
#' }

E_DHx_HmDm_HT_CIdHt.f <-
function(Hx, Hm, Dm, mHt, sHt, par.lme, Rfn=list(fn="sig2"), ...){

		if(sHt > 0){

			NHx 		= length(Hx)

			qD_u  		= rep(0,NHx)                # 	CIu E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]
			qD_o  		= rep(0,NHx)                #   CIo E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]

			pD_u        = rep(0,NHx)
			pD_o        = rep(0,NHx)

			sig2_eps    = par.lme$sig2_eps

			   CIu_DHx 	= matrix(rep(0,3*NHx),nrow=NHx,ncol=3)
			cP_CIu_DHx  = matrix(rep(0,3*NHx),nrow=NHx,ncol=3)


			for (i in 1: NHx){

				SK		= E_DHx_HmDm_HT.f(Hx = Hx[i], Hm, Dm, mHt, sHt = 0, par.lme, Rfn)
				m_DHx 	= SK$DHx; s_DHx	= sqrt(as.numeric(SK$MSE_Mean))

				qD_o[i] = m_DHx

			#	qD_1 	= m_DHx - 3*s_DHx
			#	qD_2 	= E_DHx_HmDm_HT.f( Hx = Hx[i], Hm, Dm, mHt = mw_HtT - 3.0*sd_HtT, sHt = 0, par.lme)$DHx

				qD_u[i] = 0				#	min(qD_1,qD_2)

			#   ------------------------------------------------------------------------------------------------
				pD_u[i] = Int_CdN_DHx_dHt.f(qD = qD_u[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, Rfn, nGL = 51)
			#   ------------------------------------------------------------------------------------------------
				pD_o[i]	= Int_CdN_DHx_dHt.f(qD = qD_o[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, Rfn, nGL = 51)
			#   ------------------------------------------------------------------------------------------------

				alpha = 0.025

				if(pD_u[i] > alpha){qDHx = 0}else{
				if(pD_o[i] < alpha){qDHx = qD_o[i]}else{
									qDHx = uniroot(qD.rout.f, c(qD_u[i],qD_o[i]),
														  				tol = 0.001, alpha,
																  		Hx = Hx[i], Hm, Dm, mHt, sHt,
																		par.lme = par.lme, Rfn=Rfn, nGL = 51)$root
				}}

		#		Int_CdN_DHx_dHt.f(qDHx = qD_o[i], Hx = Hx[i], Hm, Dm, mw_HtT, sd_HtT, par.lme = SK.par.lme, nGL = 51)     = Int_CdN_DHx_dHt_u
		#		pD_o[i]     = Int_CdN_DHx_dHt_o

				CIu_DHx[i,1] = m_DHx
				CIu_DHx[i,2] = qDHx
				CIu_DHx[i,3] = Hx[i]

				cP_CIu_DHx[i,1] = Int_CdN_DHx_dHt.f(qD = CIu_DHx[i,1], Hx = Hx[i], Hm, Dm,
				                                    mHt, sHt, par.lme, Rfn, nGL = 51)
				cP_CIu_DHx[i,2] = Int_CdN_DHx_dHt.f(qD = CIu_DHx[i,2], Hx = Hx[i], Hm, Dm,
				                                    mHt, sHt, par.lme, Rfn, nGL = 51)
				cP_CIu_DHx[i,3] = Hx[i]

		#		cbind(pD_u[i],pD_o[i],cP_CIu_DHx[i,1])
			}

		#	lines(Hx,CIu_DHx[,2],col=c("blue","blue","blue"), lty = c(2,1,2), lwd = c(2,2,2))

			cbind(CIu_DHx[,3],CIu_DHx[,1:2],cP_CIu_DHx[,1:2])

		#   ************************************************************************************************
		#   				Obergrenze CI(D(Hx)) fuer Schaftkurve kalibriert mit Hm,Dm und Tarifhoehe
		#   ************************************************************************************************

		#	NHx 		= 25
		#	Hx 			= seq(0,mw_HtT + 2.0*sd_HtT,length.out=NHx)    #   Hx fuer obere Grenze (SK_o:mH+2sH)

			qD_u  		= rep(0,NHx)                # 	CIu E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]
			qD_o  		= rep(0,NHx)                #   CIo E[[D(Hx)|Ht/Hm,Dm]|N(Ht|muHT(D1.3),sdHT((D1.3))]

			pD_u        = rep(0,NHx)
			pD_o        = rep(0,NHx)

			   CIo_DHx 	= matrix(rep(0,3*NHx),nrow=NHx,ncol=3)
			cP_CIo_DHx  = matrix(rep(0,3*NHx),nrow=NHx,ncol=3)

			sig2_eps    = par.lme$sig2_eps


			for (i in 1: NHx){

				SK		= E_DHx_HmDm_HT.f( Hx=Hx[i], Hm, Dm, mHt, sHt=0, par.lme=par.lme, Rfn=Rfn)
				m_DHx 	= as.numeric(SK$DHx)
				s_DHx 	= sqrt(as.numeric(SK$MSE_Mean))

				qD_u[i] = m_DHx

				qD_1 	= m_DHx + 3*s_DHx
				qD_2 	= E_DHx_HmDm_HT.f( Hx=Hx[i], Hm, Dm, mHt=mHt + 3.0*sHt, sHt=0, par.lme, Rfn=Rfn)$DHx

				qD_o[i] = max(qD_1,qD_2)

			#   Gauss-Legendre-Integration Int(-inf,+inf){P[D(Hx[i]) <= qD | Ht / Hm, Dm]dHt} :.............

			#   ------------------------------------------------------------------------------------------------
				pD_u[i] = Int_CdN_DHx_dHt.f(qD = qD_u[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, Rfn, nGL = 51)
			#   ------------------------------------------------------------------------------------------------
				pD_o[i]	= Int_CdN_DHx_dHt.f(qD = qD_o[i], Hx = Hx[i], Hm, Dm, mHt, sHt, par.lme, Rfn, nGL = 51)
			#   ------------------------------------------------------------------------------------------------

				alpha = 0.975

				if(pD_u[i] > alpha){qDHx = qD_u[i]}else{
				if(pD_o[i] < alpha){qDHx = qD_o[i]}else{
									qDHx = uniroot(qD.rout.f, c(qD_u[i] ,qD_o[i]),
														  				tol = 0.001, alpha,
																  		Hx = Hx[i], Hm, Dm, mHt, sHt,
																		par.lme, Rfn=Rfn, nGL = 51)$root
				}}

		#		Int_CdN_DHx_dHt.f(qDHx = qD_o[i], Hx = Hx[i], Hm, Dm, mw_HtT, sd_HtT, par.lme = SK.par.lme, nGL = 51)     = Int_CdN_DHx_dHt_u
		#		pD_o[i]     = Int_CdN_DHx_dHt_o

				CIo_DHx[i,1] = m_DHx
				CIo_DHx[i,2] = qDHx
				CIo_DHx[i,3] = Hx[i]

				cP_CIo_DHx[i,1] = Int_CdN_DHx_dHt.f(qD = CIo_DHx[i,1], Hx = Hx[i], Hm, Dm,
				                                    mHt, sHt, par.lme = par.lme, Rfn=Rfn, nGL = 51)
				cP_CIo_DHx[i,2] = Int_CdN_DHx_dHt.f(qD = CIo_DHx[i,2], Hx = Hx[i], Hm, Dm,
				                                    mHt, sHt, par.lme = par.lme, Rfn=Rfn, nGL = 51)
				cP_CIo_DHx[i,3] = Hx[i]
			}

			cbind(CIo_DHx[,3],CIo_DHx[,1:2],cP_CIo_DHx[,1:2])

			CI = cbind(Hx, CIu_DHx[ ,2], CIu_DHx[,1], CIo_DHx[,2], cP_CIu_DHx[,2],cP_CIo_DHx[,2])

		}else{ #  mHt=DxHx.df[Idi,"Ht"][1]gemessen

			xm = Hm/mHt
			ym = Dm

			hx = Hx[order(Hx)]/mHt
			hx = apply(cbind(1,hx),1,min)

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

			CI = cbind(Hx, CI_Mean[ ,1], CI_Mean[ ,2], CI_Mean[ ,3], rep(0.025,length(Hx)),rep(0.975,length(Hx)))

		}

	#	lines(Hx,CIu_DHx[,2],col=c("blue","blue","blue"), lty = c(2,2,2), lwd = c(2,2,2))

		colnames(CI) = c("Hx", "q_DHx_u", "DHx", "q_DHx_o", "cP_DHx_u", "cP_DHx_o")
		attr(CI, "Rfn") <- Rfn

	#   ------------------------------------------------------------------------------------------------
		return(CI)
	#   ------------------------------------------------------------------------------------------------
	}
