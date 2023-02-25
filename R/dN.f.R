#' @title evaluate Normal distribution
#' @description Internal function not usually called by users
#' @param x vector of quantiles
#' @param mw vector of means
#' @param sd vector of standard deviations
#' @param ... not currently used
#' @return numeric density of normal distribution
#' @author Edgar Kublin
#' @importFrom stats dnorm

dN.f <-
function(x, mw, sd, ...){#		pdf(N(mw,sd)) = 1/(sqrt(2*pi*sd2))*exp(-1/2*((x-mw)/sd)2)
#   ------------------------------------------------------------------------------------------------
			return(as.numeric(dnorm(x = x, mean = mw, sd = sd, log = FALSE)))
	}
