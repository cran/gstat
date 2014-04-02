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
r5to10 = rr[-unsel,]
vv = variogram(PM10~1, r5to10, width=20, cutoff = 200, tlags=0:5)
}
}
