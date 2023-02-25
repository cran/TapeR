\name{HT.par}
\alias{HT.par}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
  Height tariff parameters for estimating tree heights of unmeasured
  trees in the BWI3.
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  Height is only measured on a subset of the trees on a sample
  plots. This Height tariff is used to estimate the height of the trees
  with only a dbh measurement.
}
\usage{data(HT.par)}
\format{
  The format is:
List of 4
 $ knt.mw: num [1:16] 0 0 0 0 19.6 ...
 $ coe.mw: num [1:12] 1.3 7.28 15.1 21.75 24.39 ...
 $ knt.sd: num [1:53] 0 0 0 0 7.52 ...
 $ coe.sd: num [1:49] 0 0.618 1.376 2.142 2.486 ...
}
% \details{
% %%  ~~ If necessary, more details than the __description__ above ~~
% }
% \source{
% %%  ~~ reference to a publication or URL from which the data were
% %%  obtained ~~
%   Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem taper
%   and volume prediction method based on mixed-effects B-spline
%   regression, Eur J For Res, 132:983-997.
% }
\references{
  %%  ~~ possibly secondary sources and usages ~~
  Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem taper
  and volume prediction method based on mixed-effects B-spline
  regression, Eur J For Res, 132:983-997.  
}
\examples{
data(HT.par)
## maybe str(HT.par) ; plot(HT.par) ...
}
\keyword{datasets}
