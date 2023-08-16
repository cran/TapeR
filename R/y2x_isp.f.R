#' @title squared diameter using interpolating splines
#' @description Internal function not usually called by users
#' @param x relative height
#' @param x.grd relative heights for interpolation
#' @param y.grd diameter of taper curve at relative heights \code{x.grd} for
#' interpolation
#' @param ... not currently used
#'
#' @return squared estimated diameter based on natural interpolating spline 
#' (\code{\link[stats]{splinefun}})
#' @author Edgar Kublin
#' @importFrom stats spline

y2x_isp.f <-
function(x, x.grd, y.grd, ...){

		 y = spline(xout = c(x), x = x.grd, y = y.grd, method = "natural", ties = mean)$y
 	#	 y = predict(interpSpline( x = x.grd, y = y.grd),x=x)$y

		 return(y*y)
	}
