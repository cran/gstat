\name{plot.variogram}
\alias{plot.variogram}
\title{
Plot a Sample Variogram
}
\description{
Creates a variogram plot
}
\synopsis{
plot.variogram(x, ylim, xlim, xlab = "distance", ylab = "semivariance", multipanel = TRUE, plot.numbers = FALSE, scales, ids, ...)
}
\usage{
plot.variogram(x, ylim, xlim, xlab = "distance", ylab = "semivariance", 
    multipanel = TRUE, plot.numbers = FALSE, ids = x$id, ...)
}
\arguments{
\item{x}{ object of class "variogram", obtained from the function 
\link{variogram}, possibly containing directional or cross variograms }
\item{ylim}{ numeric vector of length 2, limits of the y-axis}
\item{xlim}{ numeric vector of length 2, limits of the x-axis}
\item{xlab}{ x-axis label }
\item{ylab}{ y-axis label }
\item{model}{ in case of a single variogram: a variogram model, as 
obtained from \link{vgm} or \link{fit.variogram}, to be drawn as 
a line in the variogram plot; in case of a set of variograms and
cross variograms: a list with variogram models }
\item{multipanel}{ logical; if TRUE, directional variograms are plotted in
different panels, if FALSE, directional variograms are plotted in the same
graph, using color, colored lines and symbols to distinguish them }
\item{plot.numbers}{ logical; if TRUE, plot number of point pairs next to
each plotted semivariance symbol }
\item{scales}{ optional argument that will be passed to \code{xyplot} in
case of the plotting of variograms and cross variograms }
\item{ids}{ ids of the data variables and variable pairs }
\item{\dots}{ any arguments that will be passed to the panel plotting functions
}
}
\value{
returns (or plots) the variogram plot
}
\references{ \url{http://www.gstat.org}}
\author{ Edzer J. Pebesma }
\seealso{
\link{variogram},
\link{fit.variogram}, 
\link{vgm}
\link{variogram.line},
}
\examples{
data(meuse)
vgm1 <- variogram(log(zinc)~1, ~x+y, meuse)
plot(vgm1)
model.1 <- fit.variogram(vgm1,vgm(1,"Sph",300,1))
plot(vgm1, model=model.1)
plot(vgm1, plot.numbers = TRUE, pch = "+")
vgm2 <- variogram(log(zinc)~1, ~x+y, meuse, alpha=c(0,45,90,135))
plot(vgm2)
# the following is only to show that directional models also work:
model.2 <- vgm(.59,"Sph",926,.06,anis=c(0,0.3))
plot(vgm2, model=model.2)
}

\keyword{dplot}

