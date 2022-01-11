#' @title construct B-Splines design matrices
#' @description Internal function not usually called by users
#' @param x relative height measurements
#' @param knt knot positions for B-Splines, usually taken from a model fit by 
#' \code{\link{TapeR_FIT_LME.f}}
#' @param ord order of B-Splines, usually taken from a model fit by 
#' \code{\link{TapeR_FIT_LME.f}}
#' @param ... not currently used
#'
#' @return List with height measurements (\code{x}), the fixed effects B-splines
#' matrix and the random effects B-splines matrix.
#' @author Edgar Kublin
#' @seealso \code{\link{TapeR_FIT_LME.f}}

XZ_BSPLINE.f <- function(x, knt, ord, ...){
    
    BS <- BSplines(knots = knt, ord = ord, der = 0, x = x)
    BS <- BS[, -ncol(BS), drop=F]
    
    return(BS)
  }
