\name{krige}
\docType{methods}
\alias{krige}
\alias{krige.locations}
\alias{krige.spatial}
\alias{idw}
\alias{idw.locations}
\alias{idw.spatial}
\alias{krige-methods}
\alias{idw-methods}
\alias{krige,formula,formula-method}
\alias{krige,formula,Spatial-method}
\alias{krige,formula,NULL-method}
\alias{idw,formula,formula-method}
\alias{idw,formula,Spatial-method}
\title{ Simple, Ordinary or Universal, global or local, Point or Block Kriging,
or simulation. }
\description{
Function for simple, ordinary or universal kriging (sometimes called
external drift kriging), kriging in a local neighbourhood, point kriging
or kriging of block mean values (rectangular or irregular blocks), and
conditional (Gaussian or indicator) simulation equivalents for all kriging
varieties, and function for inverse distance weighted interpolation. 
For multivariable prediction, see \link{gstat} and \link{predict.gstat}
}
\usage{
krige(formula, locations, ...)
krige.locations(formula, locations, data, newdata, model, ..., beta, nmax
= Inf, nmin = 0, maxdist = Inf, block, nsim = 0, indicators = FALSE,
na.action = na.pass)
krige.spatial(formula, locations, newdata, model, ..., beta, nmax
= Inf, nmin = 0, maxdist = Inf, block, nsim = 0, indicators = FALSE,
na.action = na.pass)
idw(formula, locations, ...)
idw.locations(formula, locations, data, newdata, nmax = Inf, 
	nmin = 0, maxdist = Inf, block, na.action = na.pass, idp = 2.0)
idw.spatial(formula, locations, newdata, nmax = Inf, nmin = 0, 
    maxdist = Inf, block = numeric(0), na.action = na.pass, idp = 2.0)
}
\arguments{
 \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and simple kriging use the formula \code{z~1};
  for simple kriging also define \code{beta} (see below); for universal
  kriging, suppose \code{z} is linearly dependent on \code{x} and \code{y},
  use the formula \code{z~x+y}}
  \item{locations}{ formula with only independent variables that define the
  spatial data locations (coordinates), e.g. \code{~x+y}, or object of
  class \code{Spatial}}
 \item{data}{ data frame: should contain the dependent variable, independent
  variables, and coordinates, should be missing if locations contains data. }
 \item{newdata}{ data frame or Spatial object with prediction/simulation 
  locations; should 
  contain attribute columns with the independent variables (if present) and 
  (if locations is a formula) the coordinates with names as defined in \code{locations} }
 \item{model}{ variogram model of dependent variable (or its residuals), 
  defined by a call to \link{vgm} or \link{fit.variogram}}
 \item{beta}{ only for simple kriging (and simulation based on simple
  kriging); vector with the trend coefficients (including intercept);
  if no independent variables are defined the model only contains an
  intercept and this should be the simple kriging mean }
 \item{nmax}{ for local kriging: the number of nearest observations that
  should be used for a kriging prediction or simulation, where nearest
  is defined in terms of the space of the spatial locations. By default,
  all observations are used }
 \item{nmin}{ for local kriging: if the number of nearest observations
  within distance \code{maxdist} is less than \code{nmin}, a missing 
  value will be generated; see maxdist }
 \item{maxdist}{ for local kriging: only observations within a distance
  of \code{maxdist} from the prediction location are used for prediction
  or simulation; if combined with \code{nmax}, both criteria apply }
 \item{block}{ block size; a vector with 1, 2 or 3 values containing
  the size of a rectangular in x-, y- and z-dimension respectively
  (0 if not set), or a data frame with 1, 2 or 3 columns, containing
  the points that discretize the block in the x-, y- and z-dimension
  to define irregular blocks relative to (0,0) or (0,0,0)---see also the details 
  section of \link{predict.gstat}. By default, predictions or simulations 
  refer to the support of the data values. }
 \item{nsim}{ integer; if set to a non-zero value, conditional simulation
  is used instead of kriging interpolation. For this, sequential Gaussian
  or indicator simulation is used (depending on the value of 
  \code{indicators}), following a single random path through the data.  }
 \item{indicators}{ logical, only relevant if \code{nsim} is non-zero; if
  TRUE, use indicator simulation; else use Gaussian simulation }
 \item{na.action}{ function determining what should be done with missing
  values in 'newdata'.  The default is to predict 'NA'.  Missing values 
  in coordinates and predictors are both dealt with. }
 \item{\dots}{ other arguments that will be passed to \link{gstat}}
 \item{idp}{numeric; specify the inverse distance weighting power}
}
\section{Methods}{
\describe{
\item{formula = "formula", locations = "formula"}{ 
locations specifies which coordinates in \code{data} refer to spatial coordinates
}
\item{formula = "formula", locations = "Spatial"}{ 
Object locations knows about its own spatial locations
}
\item{formula = "formula", locations = "NULL"}{ used in case of unconditional simulations;
newdata needs to be of class Spatial }
}}
\details{
Function \code{krige} is a simple wrapper method around \link{gstat}
and \link{predict.gstat} for univariate kriging prediction and conditional
simulation methods available in gstat. For multivariate prediction or
simulation, or for other interpolation methods provided by gstat (such as
inverse distance weighted interpolation or trend surface interpolation)
use the functions \link{gstat} and \link{predict.gstat} directly.

Function \code{idw} performs just as \code{krige} without a model being
passed, but allows direct specification of the inverse distance weighting
power. Don't use with predictors in the formula.

For further details, see \link{predict.gstat}.
}

\value{
a data frame containing the coordinates of \code{newdata}, and columns
of prediction and prediction variance (in case of kriging) or the
\code{abs(nsim)} columns of the conditional Gaussian or indicator
simulations }

\references{ N.A.C. Cressie, 1993, Statistics for Spatial Data,
Wiley. 

\url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers \& Geosciences, 30: 683-691.
}
\author{ Edzer J. Pebesma }
\note{  
Daniel G. Krige is a South African scientist who was a mining engineer
when he first used generalised least squares prediction with spatial
covariances in the 50's. George Matheron coined the term \code{kriging}
in the 60's for the action of doing this, although very similar approaches
had been taken in the field of meteorology. Beside being Krige's name,
I consider "krige" to be to "kriging" what "predict" is to "prediction".
}

\seealso{ \link{gstat}, \link{predict.gstat} }

\examples{
data(meuse)
data(meuse.grid)
m <- vgm(.59, "Sph", 874, .04)
# ordinary kriging:
x <- krige(log(zinc)~1, ~x+y, model = m, data = meuse, newd = meuse.grid)
library(lattice)
levelplot(var1.pred~x+y, x, aspect = "iso",
	main = "ordinary kriging predictions")
levelplot(var1.var~x+y, x, aspect = "iso",
	main = "ordinary kriging variance")
# simple kriging:
x <- krige(log(zinc)~1, ~x+y, model = m, data = meuse, newdata = meuse.grid, 
	beta=5.9)
# residual variogram:
m <- vgm(.4, "Sph", 954, .06)
# universal block kriging:
x <- krige(log(zinc)~x+y, ~x+y, model = m, data = meuse, newdata = 
	meuse.grid, block = c(40,40))
levelplot(var1.pred~x+y, x, aspect = "iso",
	main = "universal kriging predictions")
# add grid:
levelplot(var1.var~x+y, x, aspect = "iso",
		panel = function(...) {
			panel.levelplot(...)
			panel.abline(h = 0:3*1000 + 330000, v= 0:2*1000 + 179000, col = "grey")
		},
	main = "universal kriging variance")

}
\keyword{ models }
