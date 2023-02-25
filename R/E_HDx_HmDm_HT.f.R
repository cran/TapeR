#' @title Estimate height of given diameter
#' @description Calibrates a taper curve based on at least one diameter 
#' measurement and returns the height of a given diameter
#' @param Dx Scalar. Diameter for which to return height.
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
#' @details returns the height given a certain diameter.
#' @return A scalar. Estimated height (m) given a diameter.
#' @author Edgar Kublin
#' @references Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem 
#' taper and volume prediction method based on mixed-effects B-spline 
#' regression, Eur J For Res, 132:983-997.
#' @seealso \code{\link{TapeR_FIT_LME.f}}
#' @importFrom stats uniroot
#' @export
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
#' 
#' ## Calculate the height given a certain diameter threshold, say 8.5 cm
#' ht.tree1.d8.5 <- E_HDx_HmDm_HT.f (Dx=8.5, 
#'                                   Hm=tree1$Hx[3],
#'                                   Dm=tree1$Dx[3], 
#'                                   mHt = tree1$Ht[1], 
#'                                   sHt = 1, 
#'                                   par.lme = SK.par.lme,
#'                                   Rfn=list(fn="sig2")) 
#' # add to plot
#' points(x=ht.tree1.d8.5, y=8.5, pch=8, col=2, cex=2)
#' 

E_HDx_HmDm_HT.f <- function(Dx, Hm, Dm, mHt, sHt = 0, par.lme, Rfn=list(fn="sig2"), ...){

  #	Dx = D1.3
  xseq <- seq(0, mHt, length.out = 101)
  Ddiff <- E_DHx_HmDm_HT.f(xseq, Hm, Dm, mHt, sHt = 0, par.lme, Rfn)$DHx - Dx
  xnull <- xseq[which(Ddiff == 0)]
  Dprod <- Ddiff[1:100] * Ddiff[2:101]
  for(i in which(Dprod < 0)){
    xnull <- c(xnull, 
               uniroot(Hx_root.f, c(xseq[i], xseq[i + 1]), tol = 0.00001, 
                       Dx, Hm, Dm, mHt, sHt, par.lme, Rfn)$root)
  }
  
  Hx = max(xnull)
  
  return(Hx)

	}
