\name{vgmST}
\alias{vgmST}
\title{Constructing a spatio-temporal variogram}
\description{
Constructs a spatio-temporal variogram of a given type checking for a minimal set of parameters.
}
\usage{
vgmST(stModel, ..., space, time, joint, sill, nugget, stAni)
}
\arguments{
  \item{stModel}{A string indentifying the spatio-temporal variogram model.}
  \item{\dots}{unused, but ensure an exact match of the following parameters.}
  \item{space}{A spatial variogram.}
  \item{time}{A temporal variogram.}
  \item{joint}{A joint spatio-temporal variogram.}
  \item{sill}{A joint spatio-temporal sill.}
  \item{nugget}{A joint spatio-temporal nugget.}
  \item{stAni}{A spatio-temporal anisotropy.}
}

\details{
The different implemented spatio-temporal variogram models have the follwoing required parameters (see as well the example section)
\describe{
\item{separable:}{A variogram for \code{space} and \code{time} each and a joint spatio-temporal \code{sill} (variograms may have a separate nugget effect, but their joint sill will be 1) generating the call
  \preformatted{vgmST("separable", space, time, sill)}}
\item{productSum:}{A variogram (without nugget effect) for \code{space} and \code{time} each, a joint spatio-temporal \code{sill} and \code{nugget} generating the call
  \preformatted{vgmST("productSum", space, time, sill, nugget)}}
\item{sumMetric:}{A variogram (potentially including a nugget effect) for \code{space}, \code{time} and \code{joint} each and a spatio-temporal anisotropy ratio \code{stAni} generating the call
  \preformatted{vgmST("sumMetric", space, time, joint, stAni)}}
\item{simpleSumMetric:}{A variogram (without nugget effect) for \code{space}, \code{time} and \code{joint} each, a joint spatio-temporal \code{nugget} effect and a spatio-temporal anisotropy ratio \code{stAni} generating the call
  \preformatted{vgmST("simpleSumMetric", space, time, joint, nugget, stAni)}}
\item{metric:}{A spatio-temporal \code{joint} variogram (potentially inclduding a nugget effect) and \code{stAni} generating the call
  \preformatted{vgmST("metric", joint, stAni)}}}
}


\value{
Returns an S3 object  of class \code{StVariogramModel}.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{fit.StVariogram}}, \code{\link{variogramSurface}}
}
\examples{
# separable model: spatial and temporal sill will be ignored
# and kept constant at 1-nugget respectively. A joint sill is used.
separableModel <- vgmST("separable", 
                        space=vgm(0.9,"Exp", 147, 0.1),
                        time =vgm(0.9,"Exp", 3.5, 0.1),
                        sill=40)

# product sum model: spatial and temporal nugget will be ignored and kept
# constant at 0. Only a joint nugget is used.
prodSumModel <- vgmST("productSum",
                      space=vgm(39, "Sph", 343, 0),
                      time= vgm(36, "Exp",   3, 0), 
                      sill=41, nugget=17)

# sum metric model: spatial, temporal and joint nugget will be estimated
sumMetricModel <- vgmST("sumMetric",
                        space=vgm( 6.9, "Lin", 200, 3.0),
                        time =vgm(10.3, "Lin",  15, 3.6),
                        joint=vgm(37.2, "Exp",  84,11.7),
                        stAni=77.7)
                       
# simplified sumMetric model, only a overall nugget is fitted. The spatial, 
# temporal and jont nuggets are set to 0.
simpleSumMetricModel <- vgmST("simpleSumMetric",
                              space=vgm(20,"Lin", 150, 0),
                              time =vgm(20,"Lin", 10,  0),
                              joint=vgm(20,"Exp", 150, 0),
                              nugget=1, stAni=15)

# metric model
metricModel <- vgmST("metric",
                     joint=vgm(60, "Exp", 150, 10),
                     stAni=60)
}

\keyword{models}
