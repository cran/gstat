\name{plot.gstatVariogram}
\alias{plot.gstatVariogram}
\alias{plot.variogramMap}
\title{
Plot a Sample Variogram
}
\description{
Creates a variogram plot
}
\usage{
plot.gstatVariogram(x, model = NULL, ylim, xlim, xlab = "distance", 
	ylab = "semivariance", multipanel = TRUE, plot.numbers = FALSE, 
	scales, ids = x$id, group.id = TRUE, skip, ...)
plot.variogramMap(x, np = F, skip, ...)
}
\arguments{
\item{x}{ object of class "gstatVariogram", obtained from the function 
\link{variogram}, possibly containing directional or cross variograms }
\item{model}{ in case of a single variogram: a variogram model, as 
obtained from \link{vgm} or \link{fit.variogram}, to be drawn as 
a line in the variogram plot; in case of a set of variograms and
cross variograms: a list with variogram models }
\item{ylim}{ numeric vector of length 2, limits of the y-axis}
\item{xlim}{ numeric vector of length 2, limits of the x-axis}
\item{xlab}{ x-axis label }
\item{ylab}{ y-axis label }
\item{multipanel}{ logical; if TRUE, directional variograms are plotted in
different panels, if FALSE, directional variograms are plotted in the same
graph, using color, colored lines and symbols to distinguish them }
\item{plot.numbers}{ logical; if TRUE, plot number of point pairs next to
each plotted semivariance symbol }
\item{scales}{ optional argument that will be passed to \code{xyplot} in
case of the plotting of variograms and cross variograms }
\item{ids}{ ids of the data variables and variable pairs }
\item{group.id}{ logical; control for directional multivariate variograms:
if TRUE, panels divide direction and colors indicate variables
(ids), if FALSE panels divide variables/variable pairs and colors
indicate direction}
\item{skip}{ logical; can be used to arrange panels, see \code{xyplot}}
\item{np}{ logical; if TRUE, plot number of point pairs, if FALSE plot
semivariances }
\item{\dots}{ any arguments that will be passed to the panel plotting functions
(such as \code{auto.key} in examples below) }
}
\value{
returns (or plots) the variogram plot
}
\references{ \url{http://www.gstat.org}}
\author{ Edzer J. Pebesma }
\note{ currently, plotting models and/or point pair numbers is not
supported when a variogram is both directional and multivariable; also,
three-dimensional directional variograms will probably not be displayed
correctly. }
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
# the following demonstrates plotting of directional models:
model.2 <- vgm(.59,"Sph",926,.06,anis=c(0,0.3))
plot(vgm2, model=model.2)
 
g = gstat(id="zinc < 200", form=I(zinc<200)~1,loc=~x+y,data=meuse)
g = gstat(g, id="zinc < 400", form=I(zinc<400)~1,loc=~x+y,data=meuse)
g = gstat(g, id="zinc < 800", form=I(zinc<800)~1,loc=~x+y,data=meuse)
# calculate multivariable, directional variogram:
v = variogram(g, alpha=c(0,45,90,135))
plot(v, group.id = FALSE, auto.key = TRUE) # id and id pairs panels
plot(v, group.id = TRUE, auto.key = TRUE)  # direction panels

if (require(sp)) {
  plot(variogram(g, cutoff=1000, width=100, map=TRUE),
  	main = "(cross) semivariance maps")
  plot(variogram(g, cutoff=1000, width=100, map=TRUE), np=TRUE,
  	main = "number of point pairs")
}
}

\keyword{dplot}
