\name{fit.lmc}
\alias{fit.lmc}
\title{ Fit a Linear Model of Coregionalization to a Multivariable Sample Variogram }
\description{
Fit a Linear Model of Coregionalization to a Multivariable Sample Variogram;
in case of a single variogram model (i.e., no nugget) this is equivalent to
Intrinsic Correlation
}
\usage{
fit.lmc(v, g, model, fit.ranges = FALSE, fit.lmc = !fit.ranges, ...) 
}
\arguments{
\item{v}{ multivariable sample variogram, output of \link{variogram} }
\item{g}{ gstat object, output of \link{gstat} }
\item{model}{ variogram model, output of \link{vgm}; if supplied 
this value is used as initial value for each fit }
\item{fit.ranges}{ logical; determines whether the range coefficients
(excluding that of the nugget component) should be fitted; or logical
vector: determines for each range parameter of the variogram model
whether it should be fitted or fixed. }
\item{fit.lmc}{ logical; if TRUE, each coefficient matrices of partial
sills is guaranteed to be positive definite }
\item{...}{ parameters that get passed to \link{fit.variogram} }
}
\value{
returns an object of class \code{gstat}, with fitted variograms;
}
\references{ \url{http://www.gstat.org/} }
\author{ Edzer J. Pebesma }
\note{ This function does not use the iterative procedure proposed by
M. Goulard and M. Voltz (Math. Geol., 24(3): 269-286; reproduced in
Goovaerts' 1997 book) but uses simply two steps: first, each variogram
model is fitted to a direct or cross variogram; next each of the partial
sill coefficient matrices is approached by its in least squares sense
closest positive definite matrices (by setting any negative eigenvalues
to zero). }
\seealso{ \link{variogram}, \link{vgm}, \link{fit.variogram},
\code{demo(cokriging)} }
\examples{ 
}
\keyword{models}
