#' @title update variance-covariance matrix
#' @description Internal function not usually called by users
#' @param var vector of new variances 
#' @param kov variance-covariance matrix to be updated
#' @details Firstly, \code{kov} is transformed to a correlation matrix. 
#' Secondly, this correlation matrix is backtransformed to a variance-covariance
#' matrix using the new variances.
#' @return updated variance-covariance matrix
#' @author Christian Vonderach
#' @importFrom stats cov2cor

updateKOV <- function(var, kov){
  # var <- MSE_Mean
  # kov <- SK_m$KOV_Mean
  ## it might happen that variance and covariance of kov are zero
  ## suppress warning in that case
  sd <- sqrt(var)
  suppressWarnings(vcov <- outer(sd, sd) * cov2cor(kov))
  vcov[which(!is.finite(vcov))] <- 0
  
  return(vcov)
}
