\name{TapeR_FIT_LME.f}
\alias{TapeR_FIT_LME.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a taper curve model to the specified
diameter-height data.
}
\description{
Fits a taper curve model with random effects on tree-level based on B-Splines to the specified
diameter-height data. Number and position of nodes and order of B-Splines can be specified.
}
\usage{
TapeR_FIT_LME.f(Id, x, y, knt_x, ord_x, knt_z, ord_z, IdKOVb = "pdSymm",
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Id}{
Vector of tree identifiers of same length as diameter and height measurements.
}
  \item{x}{
Numeric vector of height measurements (explanatory variables) along the stem relative to
the tree height.
}
  \item{y}{
Numeric vector of diameter measurements (response) along the stem (in centimeters).
}
  \item{knt_x}{
Numeric vector of relative knot positions for fixed effects.
}
  \item{ord_x}{
Numeric scalar. Order of fixed effects Spline (4=cubic).
}
  \item{knt_z}{
Numeric vector of relative knot positions for random effects.
}
  \item{ord_z}{
Numeric scalar. Order of random effects Spline (4=cubic).
}
  \item{IdKOVb}{
Character string. Type of covariance matrix used by
\code{lme}. Only "pdSymm" makes sense. Rather reduce number of knots if
function does not converge.
}
\item{...}{Currently unused}
}
\details{
If too few trees are given, the linear mixed model (lme) will not converge. See
examples for a suggestion of node positions.

The variance parameters \code{theta} are stored in the natural parametrization (Pinheiro
  and Bates (2004), p. 93). This means log for variances and logit for
  covariances. \code{theta} is the vectorized
  triangle of the random effects covariance matrix + the residual variance
  (lSigma). Given there are 2 inner knots for random effects, the
  structure will be c(sig^2_b1, sig_b1 sig_b2, sig_b1 sig_b3, sig_b1
  sig_b4, sig^2_b2,...,sig^2_b4, lSigma)
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
List of model properties
\item{fit.lme}{Summary of the fitted lme model.}
\item{par.lme}{List of model parameters (e.g., coefficients and
  variance-covariance matrices) needed for volume estimation and other
  functions in this package.}
Components of the \code{par.lme} list
\item{knt_x}{Relative positions of the fixed effects Spline knots along the stem.}
\item{ord_x}{Order of the spline.}
\item{knt_z}{Relative positions of the random effects Spline knots along
  the stem.}
\item{ord_z}{Order of the spline.}
\item{b_fix}{Fixed-effects spline coefficients.}
\item{KOVb_fix}{Covariance of fixed-effects.}
\item{sig2_eps}{Residual variance.}%hoch2 ist problematisch da nicht ascii
\item{dfRes}{Residual degrees of freedom.}
\item{KOVb_rnd}{Covariance of random effects.}
\item{theta}{Variance parameters in natural parametrization. See Details. }
\item{KOV_theta}{Approximate asymptotic covariance matrix of variance parameters.}
%%\item{}{}
}
\references{
  Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem taper
  and volume prediction method based on mixed-effects B-spline
  regression, Eur J For Res, accepted with minor revisions (11 June
  2013).
}
\author{
Edgar Kublin
}
% \note{
% %%  ~~further notes~~
% }
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{E_DHx_HmDm_HT.f}}
}
\examples{
#load example data
data(DxHx.df)

#prepare the data (could be defined in the function directly)
Id = DxHx.df[,"Id"]
x = DxHx.df[,"Hx"]/DxHx.df[,"Ht"]#calculate relative heights
y = DxHx.df[,"Dx"]

#define the relative knot positions and order of splines
knt_x = c(0.0, 0.1, 0.75, 1.0);	ord_x = 4 # B-Spline knots: fix effects; order (cubic = 4)
knt_z = c(0.0, 0.1       ,1.0); ord_z = 4 # B-Spline knots: rnd effects

#fit the model
taper.model <- TapeR_FIT_LME.f(Id, x, y, knt_x, ord_x, knt_z, ord_z,
IdKOVb = "pdSymm")

#save model parameters for documentation or dissimination
#parameters can be load()-ed and used to predict the taper 
#or volume using one or several measured dbh
spruce.taper.pars <- taper.model$par.lme
save(spruce.taper.pars, file="spruce.taper.pars.rdata")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
