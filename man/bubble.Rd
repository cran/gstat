\name{bubble}
\alias{bubble}
\title{
Create a bubble plot of spatial data
}
\description{
Create a bubble plot of spatial data, with options for bicolour
residual plots
}
\usage{
bubble(data, xcol = 1, ycol = 2, zcol = 3, fill = TRUE, maxsize = 3,
    do.sqrt = TRUE, pch, col = c(2,3), key.entries = quantile(data[,zcol]),
    ...)
}
\arguments{
\item{data}{ data frame from which x- and y-coordinate and z-variable are taken}
\item{xcol}{ x-coordinate column number or (quoted) name }
\item{ycol}{ y-coordinate column number or (quoted) name }
\item{zcol}{ z-variable column number or (quoted) name }
\item{fill}{ logical; if TRUE, filled circles are plotted (pch = 16), 
else open circles (pch = 1); the pch argument overrides this }
\item{maxsize}{ \code{cex} value for largest circle }
\item{do.sqrt}{ logical; if TRUE the plotting symbol area (sqrt(diameter)) 
is proportional to the value of the z-variable; if FALSE, the symbol size
(diameter) is proportional to the z-variable }
\item{pch}{ plotting character }
\item{col}{ colours to be used; numeric vector of size two: first value
is for negative values, second for positive values. }
\item{key.entries}{ the values that will be plotted in the key; by default 
the five quantiles min, q.25, median q.75, max }
\item{...}{ arguments, passed to \code{xyplot}}
}
\value{
returns (or plots) the bubble plot
}
\references{ }
\author{ Edzer J. Pebesma }
\seealso{xyplot, \link{mapasp}}
\examples{
data(meuse)
bubble(meuse, max = 2.5, main = "cadmium concentrations (ppm)", 
    key.entries = c(.5,1,2,4,8,16))
bubble(meuse, "x", "y", "zinc", main = "zinc concentrations (ppm)",
    key.entries =  100 * 2^(0:4))
}
\keyword{dplot}
