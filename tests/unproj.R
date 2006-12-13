library(rgdal)
data(meuse)
coordinates(meuse) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
meuse.ll = spTransform(meuse, CRS("+proj=longlat"))
library(gstat)
variogram(log(zinc)~1, meuse.ll)

cloud1 = variogram(log(zinc)~1, meuse, cloud=T, cutoff=6000)
cloud2 = variogram(log(zinc)~1, meuse.ll, cloud=T, cutoff=6)

plot(cloud1$dist/1000, cloud2$dist, xlab="Amersfoort, km", ylab = "Long/lat")
abline(0,1)

library(fields)
data(ozone2)
oz = SpatialPointsDataFrame(ozone2$lon.lat, 
		data.frame(t(ozone2$y)), 
		proj4string=CRS("+proj=longlat"))
variogram(X870731~1,oz[!is.na(oz$X870731),])
utm16 = CRS("+proj=utm +zone=16")
oz.utm = spTransform(oz, utm16)
variogram(X870731~1,oz.utm[!is.na(oz$X870731),])
