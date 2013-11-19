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
                        space=vgm(0.8,"Exp", 140, 0.2),
                        time =vgm(0.7,"Exp",   6, 0.3),
                        sill=60)

data(vv)
fit.StVariogram(vv[1:77,], separableModel)}
}

\keyword{models}
