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
Returns a spatio-temporal variogram model, as S3 class StvariogramModel.
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
separableModel <- vgmST("separable",
                        space=vgm(0.9,"Exp", 147, 0.1),
                        time =vgm(0.9,"Exp", 3.5, 0.1),
                        sill=40)

data(vv)
fit.StVariogram(vv[1:77,], separableModel, method = "L-BFGS-B")
}

\keyword{models}
