\name{vgm}
\alias{vgm}
\alias{print.variogram.model}
\title{
Generate, or Add to Variogram Model
}
\description{
Generates a variogram model, or adds to an existing model. 
\code{print.variogram.model} prints the essence of a variogram  model.
}
\usage{
vgm(psill, model, range, nugget, add.to, anis)
print.variogram.model(x, ...)
}
\arguments{
\item{psill}{ (partial) sill of the variogram model component }
\item{model}{ model type, e.g. "Exp", "Sph", "Gau". Calling vgm() with
a missing model gives the complete list. }
\item{range}{ range of the variogram model component }
\item{nugget}{ nugget component of the variogram (this basically adds a nugget
compontent to the model) }
\item{add.to}{ a variogram model to which we want to add a component }
\item{anis}{ anisotropy parameters: }
\item{x}{ a variogram model to print }
\item{...}{ arguments that will be passed to \code{print}, e.g. \code{digits}
(see examples) }
}
\value{
an object of class variogram.model, which extends data.frame
}
\references{ }
\author{ Edzer J. Pebesma }
\seealso{
\link{print.variogram.model}, 
\link{fit.variogram}, 
\link{variogram.line}, 
link{variogram} for the sample variogram.
}
\examples{
vgm(10, "Exp", 300)
x <- vgm(10, "Exp", 300)
vgm(10, "Nug", 0)
vgm(10, "Exp", 300, 4.5)
vgm( 5, "Exp", 300, add.to = vgm(5, "Exp", 60, nugget = 2.5))
vgm(10, "Exp", 300, anis = c(30, 0.5))
vgm(10, "Exp", 300, anis = c(30, 10, 0, 0.5, 0.3))
x <- vgm(0.39527463, "Sph", 953.8942, nugget = 0.06105141)
x
print(x, digits = 3);
# to see all components, do
print.data.frame(x)
}

\keyword{models}

