\name{show.vgms}
\alias{show.vgms}
\title{
Plot Variogram Model Functions
}
\description{
Creates a trellis plot for a range of variogram models, possibly with nugget;
and optionally a set of Matern models with varying smoothness.
}
\usage{
show.vgms(min = 1e-12 * max, max = 3, n = 50, sill = 1, range = 1,
    models = as.character(vgm()[c(1:16)]), nugget = 0, kappa.range = 0.5,
	plot = TRUE)
}
\arguments{
\item{min}{ numeric; start distance value for semivariance calculation 
beyond the first point at exactly zero  }
\item{max}{ numeric; maximum distance for semivariance calculation 
and plotting }
\item{n}{ integer; number of points to calculate distance values }
\item{sill}{ numeric; (partial) sill of the variogram model }
\item{range}{ numeric; range of the variogram model }
\item{models}{ character; variogram models to be plotted }
\item{nugget}{ numeric; nugget component for variogram models }
\item{kappa.range}{ numeric; if this is a vector with more than one
element, only a range of Matern models is plotted with these kappa 
values }
\item{plot}{ logical; if TRUE, a plot is returned with the models
specified; if FALSE, the data prepared for this plot is returned }
}
\value{
returns a (Trellis) plot of the variogram models requested; see
examples. I do currently have strong doubts about the ``correctness''
of the ``Hol'' model. The ``Spl'' model does seem to need a very
large range value (larger than the study area?) to be of some value.

If plot is FALSE, a data frame with the data prepared to plot
is being returned.
}
\references{ \url{http://www.gstat.org}}
\author{ Edzer J. Pebesma }
\note{ the \code{min} argument is supplied because the variogram
function may be discontinuous at distance zero, surely when a positive
nugget is present. }
\seealso{
\link{vgm}, \link{variogram.line},
}
\examples{
show.vgms()
show.vgms(models = c("Exp", "Mat", "Gau"), nugget = 0.1)
# show a set of Matern models with different smoothness:
show.vgms(kappa.range = c(.1, .2, .5, 1, 2, 5, 10), max = 10)
}

\keyword{dplot}
