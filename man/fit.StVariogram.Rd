\name{fit.StVariogram}
\alias{fit.StVariogram}
\title{Fit a spatio-temporal sample variogram to a sample variogram}
\description{
Fits a spatio-temporal variogram of a given type to spatio-temporal sample variogram.
}
\usage{
fit.StVariogram(object, model, ..., wles = FALSE)
}
\arguments{
  \item{object}{The spatio-temporal sample variogram. Typically output from \code{\link{variogramST}}}
  \item{model}{The desired spatio-temporal model defined as a list. See Example section below.}
  \item{\dots}{arguments passed to \code{\link{optim}}}
  \item{wles}{Resiudals are weighted by the number of points in each lag class.}
}
\value{
Returns a list containing the ouput of \code{\link{optim}} and an entry stVgmFit containg the list object describing the fitted model.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{fit.variogram}}
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
                       nugget=6.7, stAni=77.7, stModel="sumMetric")

                       
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
fit.StVariogram(vv[1:77,],separableModel,method = "L-BFGS-B")$StVgmFit
}

\keyword{models}
