\name{predict.gstat}
\alias{predict.gstat}
\title{ Multivariable Geostatistical Prediction and Simulation }
\description{ The function provides the following prediction methods:
simple, ordinary, and universal kriging, simple, ordinary, and universal
cokriging, point- or block-kriging, and conditional simulation equivalents
for each of the kriging methods.  }
\usage{
predict.gstat(object, newdata, block = numeric(0), nsim = 0, 
	indicators = FALSE, BLUE = FALSE, debug.level = 1, mask, 
	na.action = na.pass, ...)
}
\arguments{
\item{object}{ object of class \code{gstat}, see \link{gstat} and
  \link{krige}}
\item{newdata}{ data frame with prediction/simulation locations; should 
  contain columns with the independent variables (if present) and the
  coordinates with names as defined in \code{locations} }
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
\item{indicators}{ logical; only relevant if \code{nsim} is non-zero; if
 TRUE, use indicator simulation, else use Gaussian simulation }
\item{BLUE}{ logical; if TRUE return the BLUE trend estimates only, 
if FALSE return the BLUP predictions (kriging) }
\item{debug.level}{ integer; set gstat internal debug level }
\item{mask}{ deprecated -- use na.action; 
logical or numerical vector; pattern with valid values
in newdata (marked as TRUE, non-zero, or non-NA); if mask is specified,
the returned data frame will have the same number and order of rows 
in newdata, and masked rows will be filled with NA's. }
\item{na.action}{ function determining what should be done with missing
values in 'newdata'.  The default is to predict 'NA'.  Missing values 
in coordinates and predictors are both dealt with. }
\item{...}{ ignored (but necessary for the S3 generic/method consistency) }
}

\details{ When a non-stationary (i.e., non-constant) mean is used, both
for simulation and prediction purposes the variogram model defined should
be that of the residual process, and not that of the raw observations.

The algorithm used by gstat for simulation random fields is the
sequential simulation algorithm. This algorithm scales well to large or
very large fields (e.g., more than $10^6$ nodes). Its power lies in using
only data and simulated values in a local neighbourhood to approximate the
conditional distribution at that location, see \code{nmax} in \link{krige}
and \link{gstat}. The larger \code{nmax}, the better the approximation,
the smaller \code{nmax}, the faster the simulation process. For selecting
the nearest \code{nmax} data or previously simulated points, gstat uses
a bucket PR quadtree neighbourhood search algorithm; see the reference
below.

For sequential Gaussian or indicator simulations, a random path through
the simulation locations is taken, which is usually done for sequential
simulations. The reason for this is that the local approximation of the
conditional distribution, using only the \code{nmax} neareast observed
(or simulated) values may cause spurious correlations when a regular
path would be followed. Following a single path through the locations,
gstat reuses the expensive results (neighbourhood selection and solution
to the kriging equations) for each of the subsequent simulations when
multiple realisations are requested.  You may expect a considerable speed
gain in simulating 1000 fields in a single call to \link{predict.gstat},
compared to 1000 calls, each for simulating a single field.

The random number generator used for generating simulations is the native
random number generator of the environment (R, S); fixing randomness by
setting seeds with \code{set.seed()} therefore works.

When mean coefficient are not supplied, they are generated as well from
their conditional distribution (MVN, BLUE estimate and covariance);
for a reference to the algorithm used see Abrahamsen and Benth,
Math. Geol. 33(6), page 742 and leave out all constraints.

Memory requirements for sequential simulation: let n be the product of
the number of variables, the number of simulation locations, and the
number of simulations required in a single call.  the gstat C function
\code{gstat_predict} requires a table of size n * 12 bytes to pass the
simulations back to R, before it can free n * 4 bytes. Hopefully, R does
not have to duplicate the remaining n * 8 bytes when the coordinates
are added as columns, and when the resulting matrix is coerced to
a \code{data.frame}.
}

\value{
a data frame containing the coordinates of \code{newdata}, and columns
of prediction and prediction variance (in case of kriging) or the
columns of the conditional Gaussian or indicator simulations 
}

\references{ 
N.A.C. Cressie, 1993, Statistics for Spatial Data, Wiley. 

\url{http://www.gstat.org/}

For bucket PR quadtrees, excellent demos are found at
\url{http://www.cs.umd.edu/~brabec/quadtree/index.html}
}
\author{ Edzer J. Pebesma }
\note{  }

\seealso{\link{gstat}, \link{krige}}

\examples{
# generate 5 conditional simulations
data(meuse)
data(meuse.grid)
v <- variogram(log(zinc)~1,~x+y, meuse)
m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
plot(v, model = m)
set.seed(131)
sim <- krige(formula = log(zinc)~1, locations = ~x+y, model = m, data 
	= meuse, newdata = meuse.grid, nmax = 15, beta = 5.9, nsim = 5)
# show all 5 simulation, using map.to.lev to rearrange sim:
levelplot(z~x+y|name, map.to.lev(sim, z=c(3:7)), aspect = mapasp(sim))

# calculate generalised least squares residuals w.r.t. constant trend:
g <- gstat(id = "log.zinc", formula = log(zinc)~1, locations = ~x+y, 
     model = m, data = meuse)
blue0 <- predict(g, newdata = meuse, BLUE = TRUE)
blue0$blue.res <- log(meuse$zinc) - blue0$log.zinc.pred
bubble(blue0, zcol = "blue.res", main = "GLS residuals w.r.t. constant")

# calculate generalised least squares residuals w.r.t. linear trend:
m <- fit.variogram(variogram(log(zinc)~sqrt(dist.m),~x+y, meuse), 
	vgm(1, "Sph", 300, 1))
g <- gstat(id = "log.zinc", formula = log(zinc)~sqrt(dist.m), 
	locations = ~x+y, model = m, data = meuse)
blue1 <- predict(g, newdata = meuse, BLUE = TRUE)
blue1$blue.res <- log(meuse$zinc) - blue1$log.zinc.pred
bubble(blue1, zcol = "blue.res", 
	main = "GLS residuals w.r.t. linear trend")

# unconditional simulation on a 100 x 100 grid
xy <- expand.grid(1:100, 1:100)
names(xy) <- c("x","y")
g.dummy <- gstat(formula = z~1, locations = ~x+y, dummy = TRUE, beta = 0,
	model = vgm(1,"Exp",15), nmax = 20)
yy <- predict(g.dummy, newdata = xy, nsim = 4)
# show one realisation:
levelplot(sim1~x+y, yy, aspect = mapasp(yy))
# show all four:
levelplot(z~x+y|name, map.to.lev(yy, z=c(3:6)), aspect = mapasp(yy))
}

\keyword{ models }
