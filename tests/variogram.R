library(gstat)
data(meuse)
variogram(log(zinc)~1, ~x+y, meuse)
