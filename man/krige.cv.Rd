\name{krige.cv}
\alias{krige.cv}
\alias{gstat.cv}
\title{ (co)kriging cross validation, n-fold or leave-one-out }
\description{ Cross validation function for simple, ordinary or universal
point kriging, kriging in a local neighbourhood.  }
\usage{
gstat.cv(object, nfold, remove.all = FALSE, verbose = FALSE, ...)
krige.cv(formula, locations, data, model = NULL, beta = NULL, nmax = Inf, 
	nmin = 0, maxdist = Inf, nfold = nrow(data), verbose = FALSE, ...)
}
\arguments{
  \item{object}{ object of class gstat; see function \link{gstat}}
  \item{nfold}{ apply n-fold cross validation; if \code{nfold} is set
  to \code{nrow(data)} (the default), leave-one-out cross validation is
  done; if set to e.g. 5, five-fold cross validation is done }
  \item{remove.all}{ logical; if TRUE, remove observations at cross validation
  locations not only for the first, but for all subsequent variables as well }
  \item{verbose}{ logical; if TRUE, progress is printed }
  \item{\dots}{ other arguments that will be passed to \link{predict.gstat}
  in case of \code{gstat.cv}, or to \link{gstat} in case of \code{krige.cv}}
  \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and simple kriging use the formula \code{z~1};
  for simple kriging also define \code{beta} (see below); for universal
  kriging, suppose \code{z} is linearly dependent on \code{x} and \code{y},
  use the formula \code{z~x+y}}
  \item{locations}{ formula with only independent variables that define the
  spatial data locations (coordinates), e.g. \code{~x+y}; if \code{data}
  is of class \code{spatial.data.frame}, this argument may be ignored, as
  it can be derived from the data }
  \item{data}{ data frame; should contain the dependent variable, independent
  variables, and coordinates. }
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
}

\details{ Leave-one-out cross validation (LOOCV) visits a data point, and
predicts the value at that location by leaving out the observed value,
and proceeds with the next data point. (The observed value is left
out because kriging would otherwise predict the value itself.) N-fold
cross validation makes a partitions the data set in N parts. For all
observation in a part, predictions are made based on the remaining N-1
parts; this is repeated for each of the N parts. N-fold cross validation
may be faster than LOOCV. }

\value{ a data frame containing the coordinates of \code{newdata}, and
columns of prediction and prediction variance of cross validated data
points, observed values, residuals, zscore (residual divided by kriging
standard error), and fold.  }

\references{ 
\url{http://www.gstat.org/}
}
\author{ Edzer J. Pebesma }
\note{ Leave-one-out cross validation seems to be much faster in plain
(stand-alone) gstat, apparently most of the effort is spent moving data 
around from R to gstat. }

\seealso{ \link{krige}, \link{gstat}, \link{predict.gstat} }

\examples{
data(meuse)
m <- vgm(.59, "Sph", 874, .04)
# five-fold cross validation:
x <- krige.cv(log(zinc)~1, ~x+y, model = m, data = meuse, nmax = 40, nfold=5)
bubble(x, z = "residual", main = "log(zinc): 5-fold CV residuals")
}
\keyword{ models }
