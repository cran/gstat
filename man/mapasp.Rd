\name{mapasp}
\alias{mapasp}
\title{ Calculate plot aspect ratio for geographic maps }
\description{ Calculate plot aspect ratio for geographic maps }
\usage{
mapasp(data, x = data$x, y = data$y)
}
\arguments{
\item{data}{ data frame }
\item{x}{ x-coordinates }
\item{y}{ y-coordinates }
}
\value{
	diff(range(y))/diff(range(x))
}
\seealso{ \link{image.data.frame}, \link{krige}}

\keyword{dplot}
