\name{krige}
\alias{krige}
\title{ Simple, Ordinary or Universal, global or local, Point or Block Kriging }
\description{
Function for simple, ordinary or universal kriging (sometimes called
external drift kriging), kriging in a local neighbourhood, point kriging
or kriging of block mean values (rectangular or irregular blocks), and
conditional (Gaussian or indicator) simulation equivalents for all kriging
varieties.
}
\usage{ krige(formula, locations, data, newdata, model, beta, nmax = Inf,
maxdist = Inf, block, nsim = 0, indicators = FALSE, ...) }
\arguments{
 \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and simple kriging use the formula \code{z~1};
  for simple kriging also define \code{beta} (see below); for universal
  kriging, suppose \code{z} is linearly dependent on \code{x} and \code{y},
  use the formula \code{z~x+y}}
 \item{locations}{ formula with only independent variables that define the
  spatial data locations (coordinates), e.g. \code{~x+y} }
 \item{data}{ data frame; should contain the dependent variable, independent
  variables, and coordinates. }
 \item{newdata}{ data frame with prediction/simulation locations; should 
  contain columns with the independent variables (if present) and the
  coordinates with names as defined in \code{locations} }
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
 \item{maxdist}{ for local kriging: only observations within a distance
  of \code{maxdist} from the prediction location are used for prediction
  or simulation; if combined with \code{nmax}, both criteria apply }
 \item{block}{ block size; a vector with 1, 2 or 3 values containing
  the size of a rectangular in x-, y- and z-dimension respectively
  (0 if not set), or a data frame with 1, 2 or 3 columns, containing
  the points that discretize the block in the x-, y- and z-dimension; 
  the latter can be used to define irregular blocks. By default,
  predictions or simulations refer to point support values. }
 \item{nsim}{ integer; if set to a non-zero value, conditional simulation
  is used instead of kriging interpolation. For this, sequential Gaussian
  or indicator simulation is used (depending on the value of 
  \code{indicators}), following a single random path through the data.  }
 \item{indicators}{ logical, only relevant if \code{nsim} is non-zero; if
  TRUE, use indicator simulation; else use Gaussian simulation }
 \item{\dots}{ other arguments that will be passed to \link{gstat}}
}
\details{
This function is a simple wrapper function around \link{gstat} and
\link{predict.gstat} for univariate kriging prediction and conditional
simulation methods available in gstat. For multivariate prediction or
simulation, or for other interpolation methods provided by gstat (such as
inverse distance weighted interpolation or trend surface interpolation)
use the functions \link{gstat} and \link{predict.gstat} directly.

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
levelplot(var1.pred~x+y, x, aspect = mapasp(x),
	main = "ordinary kriging predictions")
levelplot(var1.var~x+y, x, aspect = mapasp(x),
	main = "ordinary kriging variance")
# simple kriging:
x <- krige(log(zinc)~1, ~x+y, model = m, data = meuse, newdata = meuse.grid, 
	beta=5.9)
# residual variogram:
m <- vgm(.4, "Sph", 954, .06)
# universal block kriging:
x <- krige(log(zinc)~x+y, ~x+y, model = m, data = meuse, newdata = 
	meuse.grid, block = c(40,40))
levelplot(var1.pred~x+y, x, aspect = mapasp(x),
	main = "universal kriging predictions")
levelplot(var1.var~x+y, x, aspect = mapasp(x),
	main = "universal kriging variance")
}
\keyword{ models }
