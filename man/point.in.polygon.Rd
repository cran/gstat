\name{point.in.polygon}
\alias{point.in.polygon}
\title{ do point(s) fall in a given polygon? }
\description{
verifies for one or more points whether they fall in a given polygon
}
\usage{
point.in.polygon(point.x, point.y, pol.x, pol.y)
}
\arguments{
\item{point.x}{numerical array of x-coordinates of points}
\item{point.y}{numerical array of y-coordinates of points}
\item{pol.x}{numerical array of x-coordinates of polygon}
\item{pol.y}{numerical array of y-coordinates of polygon}
}
\value{ logical array; FALSE if a point is strictly exterior to the
polygon, TRUE if not (point is strictly interior to polygon, point is
a vertex of polygon, or point lies on the relative interior of an edge
of polygon) }
\references{ 
Uses the C function InPoly(), in gstat file polygon.c;
InPoly is Copyright 1998 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
}
\examples{
# open polygon:
point.in.polygon(1:10,1:10,c(3,4,4,3),c(3,3,4,4))
# closed polygon:
point.in.polygon(1:10,1:10,c(3,4,4,3,3),c(3,3,4,4,3))
}
\keyword{models}
