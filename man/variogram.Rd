\name{variogram}
\alias{variogram}
\alias{variogram.gstat}
\alias{variogram.formula}
\alias{variogram.default}
\alias{print.variogram}
\alias{print.variogram.cloud}
\title{
Calculate Sample or Residual Variogram or Variogram Cloud
}
\description{
Calculates the sample variogram from data, or in case of a linear model
is given, for the residuals, with options for directional, robust,
and pooled variogram, and for irregular distance intervals.  
}
\synopsis{
variogram(object, ...)
}
\usage{
variogram(object, ...)
variogram(formula, locations, data, ...)
variogram(y, locations, X, cutoff, width, alpha, beta, tol.hor, 
    tol.ver, cressie, dX, boundaries, cloud, trend.beta, debug.level, ...)
print.variogram(v, ...)
print.variogram.cloud(v, ...)
}
\arguments{
\item{object}{object of class \code{gstat}; in this form, direct
and cross (residual) variograms are calculated for all variables and
variable pairs defined in \code{object}}
\item{formula}{formula defining the response vector and (possible) 
regressors, in case of absence of regressors, use e.g. \code{z~1}}
\item{data}{data frame where the names in formula are to be found}
\item{locations}{ spatial data locations.  For variogram.formula: a
formula with only the coordinate variables in the left hand (dependent
variable) side (see examples).

For variogram.default: a matrix, with the number of rows matching
that of y, the number of columns should match the number of spatial
dimensions spanned by the data (1 (x), 2 (x,y) or 3 (x,y,z)).  }
\item{...}{any other arguments that will be passed to \link{variogram.default}}
\item{y}{ vector with responses }
\item{X}{ (optional) matrix with regressors/covariates; the number of
rows should match that of y, the number of columns equals the number
of regressors (including intercept) }
\item{cutoff}{ spatial separation distance up to which point pairs
are included in semivariance estimates }
\item{width}{ the width of subsequent distance intervals into which
data point pairs are grouped for semivariance estimates }
\item{alpha}{ direction in  plane (x,y), in positive degrees clockwise
from positive y (North): alpha=0 for direction North (increasing y),
alpha=90 for direction East (increasing x); optional a vector of
directions in (x,y) }
\item{beta}{ direction in z, in positive degrees up from the (x,y) plane; }
optional a vector of directions
\item{tol.hor}{ horizontal tolerance angle in degrees }
\item{tol.ver}{ vertical tolerance angle in degrees }
\item{cressie}{ logical; if TRUE, use Cressie's robust variogram estimate;
if FALSE use the classical method of moments variogram estimate }
\item{dX}{ include a pair of data points $y(s_1),y(s_2)$ taken at
locations $s_1$ and $s_2$ for sample variogram calculation only when
$||x(s_1)-x(s_2)|| < dX$ with and $x(s_i)$ the vector with regressors at
location $s_i$, and $||.||$ the 2-norm.  This allows pooled estimation of
within-strata variograms (use a factor variable as regressor, and dX=0.5),
or variograms of (near-)replicates in a linear model (addressing point
pairs having similar values for regressors variables) }
\item{boundaries}{ numerical vector with distance interval boundaries; 
values should be strictly increasing }
\item{cloud}{ logical; if TRUE, calculate the semivariogram cloud }
\item{trend.beta}{vector with trend coefficients, in case they are
known. By default, trend coefficients are estimated from the data.}
\item{debug.level}{ integer; set gstat internal debug level }
\item{v}{ object of class \code{variogram} or \code{variogram.cloud}
to be printed}
% \item{...}{arguments that get passed to \code{print}, e.g. digits}
}
\value{ 
an object of class "variogram" with the following fields:
\item{np}{the number of point pairs for this estimate; 
in case of a \code{variogram.cloud} see below}
\item{dist}{the average distance of all point pairs considered
for this estimate}
\item{gamma}{the actual sample variogram estimate}
\item{dir.hor}{the horizontal direction}
\item{dir.ver}{the vertical direction}
\item{id}{the combined id pair}
If cloud is TRUE: an object of class "variogram.cloud", with the field
\code{np} encoding the numbers of the point pair that contributed to a
variogram cloud estimate, as follows. The first point is found by the
integer division of np by $2^{16}$, the second point by the remainder
of that division. \link{print.variogram.cloud} shows no \code{np} field,
but does show in addition:
\item{left}{for variogram.cloud: data id (row number) of one of 
the data pair}
\item{right}{for variogram.cloud: data id (row number) of the other 
data in the pair}
}

\references{ 
Cressie, N.A.C., 1993, Statistics for Spatial Data, Wiley.

\url{http://www.gstat.org/}

}
\author{ Edzer J. Pebesma }
\note{
}
\seealso{
\link{print.variogram},
\link{plot.variogram},
\link{plot.variogram.cloud},
for variogram models: \link{vgm},
to fit a variogram model to a sample variogram: 
\link{fit.variogram}
}
\examples{
data(meuse)
# no trend:
variogram(log(zinc)~1, loc=~x+y, meuse)
# residual variogram w.r.t. a linear trend:
variogram(log(zinc)~x+y, loc=~x+y, meuse)
# directional variogram:
variogram(log(zinc)~x+y, loc=~x+y, meuse, alpha=c(0,45,90,135))
}

\keyword{models}
