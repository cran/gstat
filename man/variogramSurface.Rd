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
  A lsit specifing the spatio-temporal variogram. See the example section below.
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
# separable model: spatial and temporal sill and nugget will be ignored
# and kept constant at 1 and 0 respectively. Only a joint sill and nugget
# are used.
separableModel <- list(space=vgm(1,"Exp", 147, 0),
                       time =vgm(1,"Exp", 3.5, 0),
                       sill=40, nugget=16, stModel="separable")

# product sum model: spatial and temporal nugget will be ignored and kept
# constant at 0. Only a joint nugget is used.
prodSumModel <- list(space=vgm(39, "Sph", 343, 0),
                     time= vgm(36, "Exp",   3, 0), 
                     sill=41, nugget=17, stModel="productSum")

# sum metric model: spatial, temporal and joint nugget will be estimated
sumMetricModel <- list(space=vgm( 6.9, "Lin", 200, 3.0),
                       time =vgm(10.3, "Lin",  15, 3.6),
                       joint=vgm(37.2, "Exp",  84,11.7),
                       stAni=77.7, stModel="sumMetric")
                       
# simplified sumMetric model, only a overall nugget is fitted. The spatial, 
# temporal and jont nuggets are set to 0.
simpleSumMetricModel <- list(space=vgm(20,"Lin", 150, 0),
                             time =vgm(20,"Lin", 10,  0),
                             joint=vgm(20,"Exp", 150, 0),
                             nugget=1, stAni=15, stModel="simpleSumMetric")

# metric model
metricModel <- list(joint=vgm(60, "Exp", 150, 10),
                    stAni=60, stModel="metric")
                       
data(vv)

if(require(lattice)) {
wireframe(model~spacelag+timelag, variogramSurface(separableModel,vv))
}

# plotting of sample and model variogram
plot(vv,separableModel)

}
\keyword{models}
\keyword{spatio-temporal}
