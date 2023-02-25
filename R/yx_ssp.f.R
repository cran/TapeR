#' @title diameter using smoothing splines
#' @description Internal function not usually called by users
#' @param x relative height
#' @param x.grd relative heights for interpolation
#' @param y.grd diameter of taper curve at relative heights \code{x.grd} for
#' interpolation
#' @param ... not currently used
#'
#' @return estimated diameter based on smoothing splines 
#' (\code{\link[stats]{smooth.spline}})
#' @author Edgar Kublin
#' @importFrom stats smooth.spline predict

yx_ssp.f <-
function(x, x.grd, y.grd, ...){
  ssp = smooth.spline(x.grd, y.grd)
  return(predict(ssp, x)$y)
  }
