% $Id: krige.Rd,v 1.25 2010-01-12 12:15:48 edzer Exp $
\name{krigeST}
\docType{methods}
\alias{krigeST}
\alias{krige,formula,ST-method} % calling krigeST
% \alias{idw,formula,ST-method}

\title{ Ordinary global Spatio-Temporal Kriging}
\description{
Function for ordinary global spatio-temporal kriging on point support
}
\usage{
krigeST(formula, data, newdata, modelList, y, ..., computeVar = FALSE,
	fullCovariance = FALSE)
}
\arguments{
 \item{formula}{ formula that defines the dependent variable as a linear
  model of independent variables; suppose the dependent variable has name
  \code{z}, for ordinary and simple kriging use the formula \code{z~1};
  for simple kriging also define \code{beta} (see below); for universal
  kriging, suppose \code{z} is linearly dependent on \code{x} and \code{y},
  use the formula \code{z~x+y}}
 \item{data}{ data frame: should contain the dependent variable, independent
  variables, and coordinates, should be missing if locations contains data. }
 \item{newdata}{ data frame or Spatial object with prediction/simulation 
  locations; should 
  contain attribute columns with the independent variables (if present) and 
  (if locations is a formula) the coordinates with names as defined in 
  \code{locations} }
 \item{modelList}{ list with named elements: space, time and/or joint depending 
  on the spatio-temporal covariance family and an entry stModel. Currently 
  implemented families that may be used for stModel are \code{separable}, 
  \code{productSum}, \code{metric}, \code{sumMetric} and \code{simpleSumMetric}. 
  See the examples section in \code{\link{fit.StVariogram}} or 
  \code{\link{variogramSurface}} for details on how to define spatio-temporal 
  covariance models.}
 \item{y}{matrix; to krige multiple fields in a single step, pass data 
  as columns of matrix \code{y}. This will ignore the value of the 
  response in \code{formula}.}
 \item{\dots}{ further arguments (currently unused)}
 \item{computeVar}{logical; if TRUE, prediction variances will be returned}
 \item{fullCovariance}{logical; if FALSE a vector with prediction 
  variances will be returned, if TRUE the full covariance matrix of 
  all predictions will be returned}
}
\details{
 Function \code{krigeST} is a R implementation of the kriging function from 
 \link{gstat} using spatio-temporal covariance modells follwoing the 
 implementaion of \code{\link{krige0}}. Function \code{krigeST} offers some 
 particular methods for ordinary spatio-temporal (ST) kriging. In particular, it does not
 support block kriging or kriging in a local neighbourhood, and does not
 provide simulation.
}

\value{
 An object of the same class as \code{newdata} (deriving from 
 \code{\linkS4class{ST}}). Attributes columns contain prediction and prediction 
 variance.
}

\references{ N.A.C. Cressie, 1993, Statistics for Spatial Data,
Wiley. 

\url{http://www.gstat.org/}

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat package.
Computers \& Geosciences, 30: 683-691.
}
\author{ Edzer Pebesma, Benedikt Graeler }

\seealso{ \code{\link{krige0}}, \code{\link{gstat}}, \code{\link{predict.gstat}} }

\examples{
library(sp)
library(spacetime)
sumMetricVgm <- vgmST("sumMetric",
                      space=vgm( 4.4, "Lin", 196.6,  3),
                      time =vgm( 2.2, "Lin",   1.1,  2),
                      joint=vgm(34.6, "Exp", 136.6, 12),
                      stAni=51.7)

data(air)

rr <- rural[,"2005-06-01/2005-06-03"]
rr <- as(rr,"STSDF")

x1 <- seq(from=6,to=15,by=1)
x2 <- seq(from=48,to=55,by=1)

DE_gridded <- SpatialPoints(cbind(rep(x1,length(x2)), rep(x2,each=length(x1))), 
                            proj4string=CRS(proj4string(rr@sp)))
gridded(DE_gridded) <- TRUE
DE_pred <- STF(sp=as(DE_gridded,"SpatialPoints"), time=rr@time)
DE_kriged <- krigeST(PM10~1, data=rr, newdata=DE_pred,
                     modelList=sumMetricVgm)
gridded(DE_kriged@sp) <- TRUE
stplot(DE_kriged)
}
\keyword{ models }
