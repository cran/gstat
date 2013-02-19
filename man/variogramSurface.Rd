\name{variogramSurface}
\alias{variogramSurface}
\title{Semivariance values for a given spatio-temporal variogram model}
\description{
Generates a surface of semivariance values given a spatio-temporal variogram model (one of separable, productSum or sumMetric)
}
\usage{
variogramSurface(model, dist_grid, covariance = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A spatio-temporal variogram model generated through \code{\link{vgmST}} or \code{\link{fit.StVariogram}}.
}
  \item{dist_grid}{
  A data.frame with two columns: \code{spacelag} and \code{timelag}.
}
  \item{covariance}{
  Shall covariances instead of semivariances be returned?
}
  \item{\dots}{
Additional arguments passed on to the underlying variogram functions.
}
}
\value{
A data.frame with columns \code{spacelag}, \code{timelag} and \code{gamma}.
}
\author{
Benedikt Graeler
}
\seealso{
See \code{\link{variogramLine}} for the spatial version and \code{\link{fit.StVariogram}} for the estimation of spatio-temporal variograms.
}
\examples{
# separable model: spatial and temporal sill will be ignored
# and kept constant at 1-nugget respectively. A joint sill is used.
separableModel <- vgmST("separable", 
                        space=vgm(0.9,"Exp", 147, 0.1),
                        time =vgm(0.9,"Exp", 3.5, 0.1),
                        sill=40)

data(vv)

if(require(lattice)) {
wireframe(model~spacelag+timelag, variogramSurface(separableModel,vv))
}

# plotting of sample and model variogram
plot(vv,separableModel)

}
\keyword{models}
\keyword{spatio-temporal}
