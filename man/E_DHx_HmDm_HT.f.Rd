\name{E_DHx_HmDm_HT.f}
\alias{E_DHx_HmDm_HT.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{E_DHx_HmDm_HT_CIdHt.f}
\alias{E_HDx_HmDm_HT.f}
\alias{E_VOL_AB_HmDm_HT.f}
\title{
  %%  ~~function to do ... ~~
  Calibrates a taper curve based on at least one diameter measurement and returns
  the expected diameters and their 
  variances, volumes for parts or the complete
  stem and threshold heights given a diameter.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Calibrates a taper curve based on e.g., dbh and upper diameter
  measurements and returns the expected diameter and its approximate
  variance given relative heights
  and taper model parameters. Confidence intervals are calculated using
  exact or approximate methods. Supportive functions return the volume for parts or the complete
  stem and threshold heights given a diameter as well as variances of
  all estimates.
}
\usage{
E_DHx_HmDm_HT.f(Hx, Hm, Dm, mHt, sHt=0, par.lme, ...)
E_DHx_HmDm_HT_CIdHt.f(Hx, Hm, Dm, mHt, sHt, par.lme, ...)
E_HDx_HmDm_HT.f (Dx, Hm, Dm, mHt, sHt=0, par.lme, ...)
E_VOL_AB_HmDm_HT.f(Hm, Dm, mHt, sHt=0, A=NULL, B=NULL, iDH, par.lme, IA, nGL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hx}{
	%%     ~~Describe \code{Hx} here~~
	Numeric vector of stem heights (m) along which to return the expected diameter.
}
  \item{Hm}{
	%%     ~~Describe \code{Hm} here~~
	Numeric vector of stem heights (m) along which diameter measurements were
	taken for calibration. Can be of length 1. Must be of same length as \code{Dm}. 
}
  \item{Dm}{
	%%     ~~Describe \code{Dm} here~~
	Numeric vector of diameter measurements (cm)
	taken for calibration. Can be of length 1. Must be of same length as \code{Hm}.
}
  \item{mHt}{
	%%     ~~Describe \code{mHt} here~~
	Scalar. Tree height (m).
}
  \item{sHt}{
	%%     ~~Describe \code{sHt} here~~
	Scalar. Standard deviation of stem height. Can be 0 if height was measured
	without error.
}
  \item{par.lme}{
	%%     ~~Describe \code{par.lme} here~~
	List of taper model parameters obtained by \code{\link{TapeR_FIT_LME.f}}.
  }
    \item{Dx}{
	Scalar. Diameter for which to return height.
  }
  \item{A}{Numeric scalar defining the lower threshold of a stem section for volume
	estimation. Depends on \code{iDH}. If \code{iDH} = "D", a diameter
	(cm), if \code{iDH} = "H", a height (m). If NULL, section starts at lowest point.}
  \item{B}{Numeric scalar defining the upper threshold of a stem section for volume
	estimation. Depends on \code{iDH}. If \code{iDH} = "D", a diameter
	(cm), if \code{iDH} = "H", a height (m). If NULL, section ends at tip.}
  \item{iDH}{Character scalar. Either "D" or "H". Type of threshold for
	section volume estimation. See \code{A} or \code{B}.}
  \item{IA}{Logic scalar. If TRUE, variance calculation of height
	estimate based on 2-point distribution. If FALSE, variance calculation of height
	estimate based on Normal approximation.}
  \item{nGL}{Numeric scalar. Number of support points for numerical
	integration.} 
  \item{...}{currently unused}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  \code{E_DHx_HmDm_HT.f} returns the calibrated taper curve and
  approximate confidence intervals which can be useful for
  plotting. Uncertainty resulting from triff height estimates if tree
  height was not measured is incorporated.

  \code{E_DHx_HmDm_HT_CIdHt.f}
  returns the calibrated taper curve and "exact" confidence
  intervals. Attention: this function is somewhat time-consuming.

  \code{E_HDx_HmDm_HT.f} returns the height given a certain diameter.

  \code{E_VOL_AB_HmDm_HT.f} calculates the volume for a complete stem or sections.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
---------------------E_DHx_HmDm_HT.f---------------------

 \item{DHx}{Numeric vector of diameters (cm) (expected value) along the heights given by \code{Hx}.}
 \item{Hx}{Numeric vector of heights (m) along which to return the expected
   diameter.}
 \item{MSE_Mean}{Mean squared error for the expected value of the diameter.}
 \item{CI_Mean}{Confidence interval. Matrix of the 95\% conf. int. for the expected value of
   the diameter (cm). First column: lower limit, second column: mean, third
   column: upper limit.}
 \item{MSE_Pred}{Mean squared error for the prediction of the diameter.}
 \item{CI_Mean}{Prediction interval. Matrix of the 95\% conf. int. for the prediction of
   the diameter (cm). First column: lower limit, second column: mean, third
   column: upper limit.}

 
 ---------------------E_DHx_HmDm_HT_CIdHt.f---------------------

 Named matrix with the following columns
 
  \item{Hx}{Numeric vector of heights (m) along which to return the expected
	diameter.} 
  \item{q_DHx_u}{Lower confidence interval (cm). (95\% CI except for
	estimates close to the stem tip.)}
  \item{DHx}{Diameter estimate (cm).}
  \item{q_DHx_o}{Upper CI (cm).}
  \item{cP_DHx_u}{Probability of observations <\code{q_DHx_u}.}
  \item{cP_DHx_o}{Probability of observations <\code{q_DHx_o}.}


  ---------------------E_HDx_HmDm_HT.f---------------------

  A scalar. Estimated height (m) given a diameter.


  ---------------------E_VOL_AB_HmDm_HT.f---------------------
  \item{E_VOL}{Estimated volume (m^3).}
  \item{VAR_VOL}{Variance of the volume estimate.}
  \item{Hm}{Height of diameter measurement (m).}
  \item{Dm}{Diameter measurement (cm).}
  \item{Ht}{Tree height (m).}
  \item{Da}{Diameter at lower section threshold (cm).}
  \item{Db}{Diameter at upper section threshold (cm).}
  \item{Ha}{Height at lower section threshold (m).}
  \item{Hb}{Height at upper section threshold (m).}  
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
  Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem taper
  and volume prediction method based on mixed-effects B-spline
  regression, Eur J For Res, accepted with minor revisions (11 June
  2013).
}
\author{
  %%  ~~who you are~~
  Edgar Kublin
}
% \note{
% %%  ~~further notes~~
% }
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{TapeR_FIT_LME.f}}
}
\examples{
## Let's take a look on the taper curve of the first tree in the
## example data set if it is only calibrated using the diameter
## measurement in 2m height.

#example data
data(DxHx.df)

#taper curve parameters based on all measured trees
data(SK.par.lme)

#select data of first tree
Idi <- (DxHx.df[,"Id"] == unique(DxHx.df$Id)[1]); DxHx.df[Idi,]
tree1 <- DxHx.df[Idi,]

## Predict the taper curve based on the diameter measurement in 2 m
## height and known height 
tc.tree1 <- E_DHx_HmDm_HT.f(Hx=1:tree1$Ht[1], Hm=tree1$Hx[3],
Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 0, par.lme = SK.par.lme) 

#plot the predicted taper curve
plot(tc.tree1$Hx, tc.tree1$DHx, type="l")
#lower CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,1], lty=2)
#upper CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,3], lty=2)
#lower prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,1], lty=3)
#upper prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,3], lty=3)
#add measured diameter
points(tree1$Hx[3], tree1$Dx[3], pch=3, col=2)
#add the observations
points(tree1$Hx, tree1$Dx)

## Add the population average taper curve (without calibration) to the
## plot (not of high practical interest but good to know how to get it).
Ht  = tree1$Ht[1]
Hx  = tree1$Hx
#get fixed-effects design matrix for the Splines
X = TapeR:::XZ_BSPLINE.f(x=Hx/Ht, par.lme = SK.par.lme)$X
#Calculate population average taper curve
DHx_PA = X \%*\% SK.par.lme$b_fix
#add to plot
lines(tree1$Hx, DHx_PA, lwd=2, lty=4)

## Let's see how CI's change if there's some uncertainty in the height
## measurement 
tc.tree1 <- E_DHx_HmDm_HT.f(Hx=1:tree1$Ht[1], Hm=tree1$Hx[3],
Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 1, par.lme = SK.par.lme) 

#plot the predicted taper curve
plot(tc.tree1$Hx, tc.tree1$DHx, type="l", xlab="Height (m)",
ylab="Diameter (cm)")
#lower CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,1], lty=2)
#upper CI
lines(tc.tree1$Hx, tc.tree1$CI_Mean[,3], lty=2)
#lower prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,1], lty=3)
#upper prediction interval
lines(tc.tree1$Hx, tc.tree1$CI_Pred[,3], lty=3)
#add measured diameter
points(tree1$Hx[3], tree1$Dx[3], pch=3, col=2)
#add the observations
points(tree1$Hx, tree1$Dx)

\dontrun{
## Calculate "exact" CIs. Careful: This takes a while!
library(pracma)# for numerical integration with gaussLegendre()
tc.tree1.exact <- E_DHx_HmDm_HT_CIdHt.f(Hx=1:tree1$Ht[1],
Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 1, par.lme =
SK.par.lme) 
#add exact confidence intervals to approximate intervals above - fits
#quite well
lines(tc.tree1.exact[,1], tc.tree1.exact[,2], lty=2,col=2)
lines(tc.tree1.exact[,1], tc.tree1.exact[,4], lty=2,col=2)


## Calculate the height given a certain diameter threshold, say 8.5 cm
ht.tree1.d8.5 <- E_HDx_HmDm_HT.f (Dx=8.5, Hm=tree1$Hx[3],
Dm=tree1$Dx[3], mHt = tree1$Ht[1], sHt = 1, par.lme = SK.par.lme) 
#add to above created plot
abline(v=ht.tree1.d8.5)

## Calculate the timber volume
#for the whole stem
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme)

#Calculate the timber volume for a selected section given a height
#threshold (0.3 - 5 m)
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, A=0.3, B=5, iDH = "H")

#Calculate the timber volume for a selected section given a diameter
#threshold (30 cm - 15 cm) (negative value if A<B)
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, A=30, B=15, iDH = "D")

#The variance estimate resulting from the tree height uncertainty using
#a Normal approximation takes much longer...
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, IA=FALSE)
proc.time() - ptm


#... than the calculation using a 2-point distribution...
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 1, par.lme = SK.par.lme, IA=TRUE)
proc.time() - ptm

#...fastest if no height variance is assumed
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 0, par.lme = SK.par.lme, IA=FALSE)
proc.time() - ptm

#Also the number of supportive points for the numerical integration
#influences the calculation time
ptm <- proc.time()
E_VOL_AB_HmDm_HT.f(Hm=tree1$Hx[3], Dm=tree1$Dx[3], mHt = tree1$Ht[1]
, sHt = 0, par.lme = SK.par.lme, IA=FALSE, nGL=10)
proc.time() - ptm
}
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
