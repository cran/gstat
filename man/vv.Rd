\name{vv}
\alias{vv}
\title{Precomputed variogram for PM10 in data set air}
\description{
Precomputed variogram for PM10 in data set air
}
\format{
data set structure is explained in \link{variogramST}.
}
\usage{
data(vv)
}
\examples{
\dontrun{
# obtained by:
data(air)
rr = rural[,"2005::2010"]
unsel = which(apply(as(rr, "xts"), 2, function(x) all(is.na(x))))
rrr = rr[-unsel,]
vv = variogram(PM10~1, rrr, width=20, cutoff = 200)
}
}
