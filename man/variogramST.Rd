\name{variogramST}
\alias{variogramST}
\title{Calculate Spatio-Temporal Sample Variogram}
\description{
Calculates the sample variogram from spatio-temporal data.
}
\usage{
variogramST(formula, locations, data, ..., tlags = 0:15,  cutoff,
	width = cutoff/15, boundaries = seq(0, cutoff, width),
	progress = interactive(), pseudo = TRUE, assumeRegular=FALSE, na.omit=FALSE)
}
\arguments{
\item{formula}{formula, specifying the dependent variable.}
\item{locations}{A STFDF or STSDF containing the variable; kept for 
  compatibility reasons with variogram, either \code{locations} or \code{data} 
  must be provided.}
\item{data}{A STFDF or STSDF containing the variable.}
\item{...}{any other arguments that will be passed to the underlying 
  \code{\link{variogram}} function.}
\item{tlags}{ integer; time lags to consider; see details below.}
\item{cutoff}{ spatial separation distance up to which point pairs
are included in semivariance estimates; as a default, the
length of the diagonal of the box spanning the data is divided by three. }
\item{width}{ the width of subsequent distance intervals into which
data point pairs are grouped for semivariance estimates, by default the 
  \code{cutoff} is divided into 15 equal lags.}
\item{boundaries}{ numerical vector with distance interval upper boundaries; 
values should be strictly increasing }
\item{progress}{ logical; if TRUE, show text progress bar }
\item{pseudo}{ integer; use pseudo cross variogram for computing
time-lagged spatial variograms? -1: find out from coordinates -- if they
are equal then yes, else no; 0: no; 1: yes. }
\item{assumeRegular}{logical; whether the time series should be assumed regular. 
  The first time step is assumed to be representative for the whole series. Note,
  that temporal lags are considered by index, and no check is made whether pairs
  actually have the desired separating distance.}
\item{na.omit}{shall all \code{NA} values in the spatio-temporal variogram be dropped? In case where complete rows or columns in the variogram consists of \code{NA} only, \code{plot} might produce a distorted picture.}
}

\value{The spatio-temporal sample variogram contains besides the fields 
  \code{np}, \code{dist} and \code{gamma} the spatio-temporal fields, 
  \code{timelag} and \code{spacelag}, the first of which indicates the time lag 
  used, the second the space lag. \code{spacelag} is the midpoint in the spatial
  lag intervals, whereas \code{dist} is the average distance between the point
  pairs found in a distance interval. To compute variograms for space lag $h$ 
  and time lag $t$, the pseudo cross variogram $(Z_i(s)-Z_{i+t}(s+h))^2$ is 
  averaged over all time lagged observation sets $Z_i$ and $Z_{i+t}$ available
  (weighted by the number of pairs involved).
}

\references{ 
Cressie, N.A.C., 1993, Statistics for Spatial Data, Wiley.

Cressie, N., C. Wikle, 2011, Statistics for Spatio-temporal Data, Wiley.

\url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers \& Geosciences, 30: 683-691.
}
\author{ Edzer Pebesma }

\seealso{
\code{\link{plot.StVariogram}},
for variogram models: \code{\link{vgmST}},
to fit a spatio-temporal variogram model to a spatio-temporal sample variogram: 
\code{\link{fit.StVariogram}}
}
\examples{
# The following spatio-temporal variogram has been calcualted through
# vv = variogramST(PM10~1, r5to10, width=20, cutoff = 200)
# in the vignette "st".

data(vv)
str(vv)
plot(vv)
}

\keyword{models}
