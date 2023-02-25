#' @title estimate cubic smoothing spline with linear extrapolation
#' @description Internal function not usually called by users
#' @param knt knots position of B-Splines
#' @param coe estimated coefficient for B-Splines
#' @param x position at which to evaluate B-Splines model
#' @param ... not currently used
#'
#' @return expected diameter given knots and coefficients at position (height)
#' \code{x}.
#' @author Edgar Kublin
#' @import splines

EYx_ssp.f <-
function(knt,coe,x, ...){
#   ************************************************************************************************

#   Berechnung kubischer smoothing.spline mit linearer Extrapolation ausserhalb der Knoten

		#	knt = knt_ssp.mw; x = min(x.plt); coe = coe_ssp.mw; x = -1

		if(x<min(knt)){# lineare Extrapolation unterhalb des kleinsten Knoten
			x_u 	= min(knt)
			X 		= splineDesign(knots = knt, x = x_u, ord = 4, derivs = c(0), outer.ok = T)
			y_u		= X%*%coe

			X  		= splineDesign(knots = knt, x = x_u, ord = 4, derivs = c(1), outer.ok = T)
			dy_u    = X%*%coe

			y_x 	= y_u + (x-x_u)*dy_u

		}else{
			 if (x>max(knt)){# lineare Extrapolation oberhalb des groessten Knoten

				x_o 	= max(knt)

				X 		= splineDesign(knots = knt, x = x_o, ord = 4, derivs = c(0), outer.ok = T)
				y_o		= X%*%coe

				X  		= splineDesign(knots = knt, x = x_o, ord = 4, derivs = c(1), outer.ok = T)
				dy_o    = X%*%coe

				y_x 	= y_o + (x-x_o)*dy_o
			}else{
				X  		= splineDesign(knots = knt, x = x, ord = 4, derivs = c(0), outer.ok = T)
				y_x    = X%*%coe
			}
		}

		return(as.numeric(y_x))
	}
