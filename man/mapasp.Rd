\name{mapasp}
\alias{mapasp}
\title{ Calculate plot aspect ratio for geographic maps }
\description{ Calculate plot aspect ratio for geographic maps }
\usage{
mapasp(data, x = data$x, y = data$y, get.number = FALSE)
}
\arguments{
\item{data}{ data frame }
\item{x}{ x-coordinates }
\item{y}{ y-coordinates }
\item{get.number}{ logical; if true the aspect value is returned else,
and in R >= 2.0.0, "iso" is returned }
}
\value{
	on R version 2.0.0 or higher, returns "iso", else returns
	diff(range(y))/diff(range(x))
}
\seealso{ \link{image.data.frame}, \link{krige}}

\keyword{dplot}
