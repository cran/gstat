\name{bpy.colors}
\alias{bpy.colors}
\title{ blue-pink-yellow color scheme }
\description{ Create a vector of `n' ``contiguous'' colors.  }
\usage{ bpy.colors(n) }
\arguments{
\item{n}{number of colors (>= 1) to be in the palette}
}
\value{ A character vector, `cv', of color names.  This can be used
either to create a user-defined color palette for subsequent graphics
by `palette(cv)', a `col=' specification in graphics functions or in
`par'.  }
\references{ see url; gnuplot has this color map }
\author{ unknown }
\url{ http://www.ihe.uni-karlsruhe.de/mitarbeiter/vonhagen/palette.en.html }
\note{ in contrast to e.g. rainbow, this color map prints well on 
black-and-white printers. }
\seealso{\link{rainbow}, \link{cm.colors}}
\examples{
bpy.colors(10)
p <- expand.grid(x=1:30,y=1:30)
p$z <- p$x + p$y
image(p, col = bpy.colors(100))
}
\keyword{color}
