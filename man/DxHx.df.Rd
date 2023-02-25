\name{DxHx.df}
\alias{DxHx.df}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
  Example dataset of 10 trees with 10 diameter and height measurements
  for each tree.
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  Example dataset of 10 trees with 10 diameter and height measurements
  for each tree.
}
\usage{data(DxHx.df)}
\format{
  A data frame with 172 observations on the following 4 variables.
  \describe{
    \item{\code{Id}}{Numeric vector of tree IDs.}
    \item{\code{Dx}}{Numeric vector of diameter measurements.}
    \item{\code{Hx}}{Numeric vector of height measurements.}
    \item{\code{Ht}}{Numeric vector of tree height (repeated for each
	  measurement in each tree).}
  }
}
\details{
  %%  ~~ If necessary, more details than the __description__ above ~~
  Measured for BWI3.
}
% \source{
% %%  ~~ reference to a publication or URL from which the data were
% %%  obtained ~~
% 
% }
\references{
  Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem taper
  and volume prediction method based on mixed-effects B-spline
  regression, Eur J For Res, 132:983-997.
}
\examples{
data(DxHx.df)
head(DxHx.df)
}
\keyword{datasets}
