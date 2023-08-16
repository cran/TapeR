#' @title builds B-Splines Matrix with appropriate knots for taper fitting
#' @description Internal function not usually called by users
#' @param knots knot positions for spline function
#' @param ord order of the spline function
#' @param der derivatives
#' @param x height measurements
#' @param ... not currently used
#' @details internally \code{\link[splines]{splineDesign}} is called
#' @return B-Splines matrix build using \code{\link[splines]{splineDesign}}
#' @author Edgar Kublin
#' @usage BSplines(knots = c(seq(0, 1, 0.1)), ord = 4, der = 0, x = c(seq(0, 1, 0.01)), ...)
#' @import splines

BSplines <- function(knots=c(seq(0, 1, 0.1)), ord=4, der=0, x=c(seq(0, 1, 0.01)), ...){

  # nur abhaengig von TapeR-Modell und der relativen Messhoehe 'x'

	TK = TransKnots(knots=knots, ord=ord) 
	BSplines = splineDesign(knots=TK, x=x, ord=ord, derivs=c(rep(der,length(x))),
	                        outer.ok=TRUE)
  return(BSplines)
}

