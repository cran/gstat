\name{plot.variogram.cloud}
\alias{plot.variogram.cloud}
\title{
Plot and Identify Data Pairs on Sample Variogram Cloud
}
\description{
Plot a sample variogram cloud, possibly with identification of
individual point pairs
}
\usage{
plot.variogram.cloud(x, identify = FALSE, xlim, ylim, xlab, ylab, ...)
}
\arguments{
\item{x}{ object of class variogram.cloud }
\item{identify}{ logical; if TRUE, the plot allows identification of
the indexes of point pairs that correspond to individual variogram
cloud points }
\item{xlim}{ limits of x-axis }
\item{ylim}{ limits of y-axis }
\item{xlab}{ x axis label }
\item{ylab}{ y axis label }
\item{...}{ parameters that are passed through to \link{plot.variogram}
(in case of identify = FALSE) or to plot (in case of identify = TRUE) }
}
\value{
if identify is T, a data frame of class \code{point.pairs} with in its 
rows the point pairs identified, if identify is F, a plot of the 
variogram cloud (see \link{plot.variogram})
}
\references{ 
\url{http://www.gstat.org/}
}
\author{ Edzer J. Pebesma }
\seealso{
\link{plot.point.pairs},
\link{variogram.formula},
\link{variogram},
\link{plot.variogram},
}
\examples{
data(meuse)
# no trend:
plot(variogram(log(zinc)~1, loc=~x+y, data=meuse, cloud=TRUE))
# plot(variogram(log(zinc)~1, loc=~x+y, data=meuse, cloud=TRUE, cressie=TRUE))
}

\keyword{dplot}
