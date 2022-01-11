#' @title Functions to put different size of uncertainty on given measurements
#' @description When estimating a tree specific taper curve based on given
#' measurements one can modify the assumed measurement uncertainty by these
#' functions
#' @param relH relative tree height for which the assumed residual variance
#' should be calculated
#' @param fn name of function to be applied as character string
#' @param sig2 residual variance from fitted model, cf.
#' \code{\link{TapeR_FIT_LME.f}}
#' @param par either NULL or a list with parameters to the different functions.
#' See details.
#' @details When estimating a tree specific taper curve based on given
#' measurements the residual variance of the model is taken into account to
#' estimate the tree specific random effects. Alternatively, it is possible to
#' make assumptions about the measurement error, eventually at specific
#' relative heights. With that, one can e.g. force the taper curve through the
#' given measurements. Standard behaviour not necessarily leads to passing the
#' measurements, if more than one is given.
#'
#' Different functions are available. \code{sig2} applies the model residual
#' variance and hence is the default behaviour.
#' \code{zero} means assuming no residual variance and forcing the taper curve
#' through the given measurements. Care has to be taken in this case because
#' forcing the taper curve through a lot of measurements might result in
#' implausible results.
#' \code{linear} interpolates between zero and the given residual variance along
#' the stem, i.e. from bottom to tree top.
#' \code{bilinear} puts zero variance not at zero but at a predefined location
#' (can be given via \code{par}). Below and above a linear interpolation is
#' done up to the given residual variance. If zero variance position is not
#' given, it is set at 5\% of tree height (approximately height of dbh).
#' \code{laglinear} assumes zero variance up to a predefined location (defaults
#' to 5\% of tree height) and interpolates upwards to the given residual
#' variance of the model.
#' \code{quadratic} function distributes residual variance according to a
#' quadratic function along the stem. It is build so that zero variance is put
#' at a predefined location (defaults to 5\% of tree height) and model residual
#' variance (as a default) at tree top.
#' \code{dnorm} and \code{dlnorm} put residual variance in form of an inverse
#' normal or an inverse log-normal distribution along the stem with a
#' zero-minimum at a predefined location (defaults to 5\% of tree height).
#' See examples for a visualisation.
#'
#' For all functions (except \code{zero} and \code{sig2}) the point of zero
#' residual variance is defined by \code{par$zrv} if given, otherwise set to 0.05.
#' For \code{dnorm} one can additionally provide the parameter \code{sd} to
#' determine standard deviation. By default it is set to \code{zrv/3}; in case
#' of \code{dlnorm} one can define \code{lsd} (sdlog, cf. \code{\link{dlnorm}}),
#' which is by default set to \code{1-sqrt(zrv)}. It is up to the user to define
#' meaningful parameters and use the functions in appropriate context.
#'
#' @return vector of assumed residual variance
#' @importFrom stats approx dnorm dlnorm
#' @author Christian Vonderach
#' @export
#'
#' @examples
#' curve(resVar(relH=x, fn = "sig2", 0.5))
#' curve(resVar(relH=x, fn = "zero", 0.5))
#' curve(resVar(relH=x, fn = "linear", 0.5))
#' curve(resVar(relH=x, fn = "bilinear", 0.5))
#' curve(resVar(relH=x, fn = "laglinear", 0.5))
#' curve(resVar(relH=x, fn = "quadratic", 0.5))
#' curve(resVar(relH=x, fn = "dnorm", 0.5))
#' curve(resVar(relH=x, fn = "dnorm", 0.5, par=list(zrv=0.2, sd=0.2/3)))
#' curve(resVar(relH=x, fn = "dlnorm", 0.5))
#' curve(resVar(relH=x, fn = "dlnorm", 0.5, par=list(zrv=0.2)))
#' invisible(sapply(seq(0.01, 0.99, length.out=20), function(a){
#'   curve(resVar(relH=x, fn = "dlnorm", 0.5, par=list(zrv=a, lsd=(1-sqrt(a)))),
#'   n=1000)
#' }))



resVar <- function(relH, fn, sig2, par=NULL){
  stopifnot(relH>=0 & relH <=1)
  stopifnot(fn %in% c("sig2", "zero", "linear", "bilinear", "laglinear",
                      "quadratic", "dnorm", "dlnorm"))
  stopifnot(sig2 > 0)
  if(is.null(par)){
    zrv <- 0.05
  } else {
    zrv <- par$zrv[1] # location of zero residual variance
  }

  if(fn[1] == "sig2"){

    return(rep(sig2, length(relH)))

  } else if(fn[1] == "zero"){

    return(rep(0, length(relH)))

  } else if(fn[1] == "linear"){

    rv <- approx(x = c(0, 1), y=c(0, sig2), xout = relH)$y
    return(rv)

  } else if(fn[1] == "bilinear"){

    rv <- ifelse(relH < zrv,
                 approx(x = c(0, zrv), y=c(sig2, 0), xout = relH)$y,
                 approx(x = c(zrv, 1), y=c(0, sig2), xout = relH)$y)
    return(rv)

  } else if(fn[1] == "laglinear"){

    rv <- approx(x = c(zrv, 1), y=c(0, sig2), xout = relH, rule = c(2, 1))$y
    return(rv)

  } else if(fn[1] == "quadratic"){

    a <- sig2 / (1 - zrv)^2
    rv <- a * (relH - zrv)^2
    return(rv)

  } else if(fn[1] == "dnorm"){

    if(is.null(par$sd)){
      sd <- zrv/3
    } else {
      sd <- par$sd
    }
    scaling <- dnorm(zrv, zrv, sd)
    rv <- (scaling - dnorm(relH, zrv, sd)) / scaling * sig2
    return(rv)

  } else if(fn[1] == "dlnorm"){

    if(is.null(par$lsd)){
      lsd <- 1 - sqrt(zrv)
    } else {
      lsd <- par$lsd
    }
    scaling <- dlnorm(zrv, log(zrv) + lsd^2, lsd)
    rv <- (scaling - (dlnorm(relH, log(zrv) + lsd^2, lsd))) / scaling * sig2
    return(rv)

  } else {

    return(rep(NA, length(relH)))

  }
}
