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
  \item{model}{The desired spatio-temporal model defined through \code{\link{vgmST}}.}
  \item{\dots}{arguments passed to \code{\link{optim}}}
  \item{wles}{Resiudals are weighted by the number of points in each lag class.}
}
\value{
Returns a spatio-temporal variogram model, as S3 class StVariogramModel. It carries the temporal and spatial unit as attributes in order to allow \code{\link{krigeST}} to adjust for different units. The units are obtained from the provided empirical variogram.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{fit.variogram}}
}
\examples{
# separable model: spatial and temporal sill will be ignored
# and kept constant at 1-nugget respectively. A joint sill is used.
\dontrun{
separableModel <- vgmST("separable",
                        space=vgm(0.9,"Exp", 123, 0.1),
                        time =vgm(0.9,"Exp", 2.9, 0.1),
                        sill=100)

data(vv)
separableModel <- fit.StVariogram(vv, separableModel,
                                  method="L-BFGS-B",
                                  lower=c(10,0,0.01,0,1),
                                  upper=c(500,1,20,1,200))
plot(vv, separableModel)
}
}

\keyword{models}
