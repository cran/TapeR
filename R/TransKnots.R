#' @title transform knot vector
#' @description Internal function not usually called by users
#' @param knots knot positions for spline function
#' @param ord order of the spline function
#' @param ... not currently used
#'
#' @return transformed knots vector, especially with repeated first and last 
#' knot given order of spline function
#' @author Edgar Kublin

TransKnots <-
function(knots=c(seq(0,1,0.1)),ord=4, ...){
#   ***********************************************************************************************
      c(rep(min(knots),ord),knots[2:(length(knots)-1)],rep(max(knots),ord))
    }
