\name{gstat}
\alias{gstat}
\alias{print.gstat}
\title{ Creates gstat Objects }
\description{
Function that creates gstat objects; objects that hold all the information
necessary for univariate or multivariate geostatistical prediction
(simple, ordinary or universal (co)kriging), or its conditional or
unconditional Gaussian or indicator simulation equivalents.  }
\usage{
gstat(g, id, formula, locations, data, model = NULL, beta, nmax = Inf,
	maxdist = Inf, dummy = FALSE, set, fill.all = FALSE, 
	variance = "identity")
print.gstat(x, ...)
}
\arguments{
  \item{g}{ gstat object to append to; if missing, a new gstat object
  is created }
  \item{id}{ id of new variable; if missing, \code{varn} is used with
    \code{n} the number for this variable. If a cross variogram is entered,
    \code{id} is a vector with the two \code{id} values , e.g. 
    \code{c("zn", "cd")} and further only supply arguments \code{g} 
    and \code{model}}
  \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and simple kriging use the formula \code{z~1};
  for simple kriging also define \code{beta} (see below); for universal
  kriging, suppose \code{z} is linearly dependent on \code{x} and \code{y},
  use the formula \code{z~x+y}}
  \item{locations}{ formula with only independent variables that define the
  spatial data locations (coordinates), e.g. \code{~x+y} }
  \item{data}{ data frame; contains the dependent variable, independent
  variables, and locations. }
  \item{model}{ variogram model for this \code{id}; defined by a call to 
  \link{vgm}; see argument \code{id} to see how cross variograms are entered }
  \item{beta}{ only for simple kriging (and simulation based on simple
  kriging); vector with the trend coefficients (including intercept);
  if no independent variables are defined the model only contains an
  intercept and this should be the simple kriging mean }
  \item{nmax}{ for local kriging: the number of nearest observations that
  should be used for a kriging prediction or simulation, where nearest
  is defined in terms of the space of the spatial locations }
  \item{maxdist}{ for local kriging: only observations within a distance
  of \code{maxdist} from the prediction location are used for prediction
  or simulation; if combined with \code{nmax}, both criteria apply }
  \item{dummy}{ logical; if TRUE, consider this data as a dummy variable
  (only necessary for unconditional simulation) }
  \item{set}{ named list with optional parameters to be passed to
  gstat (only \code{set} commands of gstat are allowed; see gstat manual) }
  \item{x}{ gstat object to print }
  \item{fill.all}{ logical; if TRUE, fill all of the variogram and cross
  variogram model slots in \code{g} with the given variogram model }
  \item{variance}{ character; variance function to transform to non-stationary
  covariances; "identity" does not transform, other options are "mu" (poisson)
  and "mu(1-mu)" (binomial) }
  \item{...}{ arguments that are passed to the printing of the variogram
  models only}
}
\details{ to print the full contents of the object \code{g} returned,
use \code{as.list(g)} }
\value{
an object of class \code{gstat}, which inherits from \code{list}.
Its components are: 
 \item{data}{list; each element is a list with the \code{formula}, 
  \code{locations}, \code{data}, \code{nvars}, and \code{beta} for a 
  variable}
 \item{model}{list; each element contains a variogram model; names are
  those of the elements of \code{data}; cross variograms have names of
  the pairs of data elements, separated by a \code{.} (e.g.: 
  \code{var1.var2}}
 \item{set}{list; named list, corresponding to set \code{name}=\code{value};
 gstat commands (look up the set command in the gstat manual for a full list)}
}
\references{ \url{http://www.gstat.org/} }
\author{ Edzer J. Pebesma }
\note{ 
the function currently copies the data objects into the gstat object, so
this may become a large object. I would like to copy only the name of the
data frame, but could not get this to work. Any help is appreciated.  
}


\seealso{ \link{predict.gstat}, \link{krige} }

\examples{
data(meuse)
# let's do some manual fitting of two direct variograms and a cross variogram
g <- gstat(id = "ln.zinc", formula = log(zinc)~1, locations = ~x+y, 
	data = meuse)
g <- gstat(g, id = "ln.lead", formula = log(lead)~1, locations = ~x+y, 
	data = meuse)
# examine variograms and cross variogram:
plot(variogram(g))
# enter direct variograms:
g <- gstat(g, id = "ln.zinc", model = vgm(.55, "Sph", 900, .05))
g <- gstat(g, id = "ln.lead", model = vgm(.55, "Sph", 900, .05))
# enter cross variogram:
g <- gstat(g, id = c("ln.zinc", "ln.lead"), model = vgm(.47, "Sph", 900, .03))
# examine fit:
plot(variogram(g), model = g$model, main = "models fitted by eye")
# see also demo(cokriging) for a more efficient approach

# Inverse distance interpolation with inverse distance power set to .5:
# (kriging variants need a variogram model to be specified)
data(meuse)
data(meuse.grid)
meuse.gstat <- gstat(id = "zinc", formula = zinc ~ 1, locations = ~ x + y,
	data = meuse, nmax = 7, set = list(idp = .5))
meuse.gstat
z <- predict(meuse.gstat, meuse.grid)
levelplot(zinc.pred~x+y, z, aspect = mapasp(z))
# see demo(cokriging) and demo(examples) for further examples, 
# and the manuals for predict.gstat and image
}
\keyword{ models }
