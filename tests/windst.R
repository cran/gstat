suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(spacetime))
suppressPackageStartupMessages(library(gstat))
suppressPackageStartupMessages(library(stars))

data(wind)
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y
proj4string(wind.loc) = "+proj=longlat +datum=WGS84 +ellps=WGS84"

wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
wind$jday = as.numeric(format(wind$time, '%j'))
stations = 4:15
windsqrt = sqrt(0.5148 * wind[stations]) # knots -> m/s
Jday = 1:366
daymeans = colMeans(
	sapply(split(windsqrt - colMeans(windsqrt), wind$jday), colMeans))
meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
velocities = apply(windsqrt, 2, function(x) { x - meanwind })
# match order of columns in wind to Code in wind.loc;
# convert to utm zone 29, to be able to do interpolation in
# proper Euclidian (projected) space:
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])
pts = SpatialPoints(pts)
#if (require(rgdal) && require(maps)) {
proj4string(pts) = "+proj=longlat +datum=WGS84 +ellps=WGS84"
utm29 = CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84")
pts = spTransform(pts, utm29)
# note the t() in:
w = STFDF(pts, wind$time, data.frame(values = as.vector(t(velocities))))

library(mapdata)
library(maptools)
m = map2SpatialLines(
    map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84 +ellps=WGS84"
m = spTransform(m, utm29)

# setup grid
grd = SpatialPixels(SpatialPoints(makegrid(m, n = 300)),
	proj4string = proj4string(m))
# grd$t = rep(1, nrow(grd))
#coordinates(grd) = ~x1+x2
#gridded(grd)=TRUE

# select april 1961:
w = w[, "1961-04"]

covfn = function(x, y = x) { 
	du = spDists(coordinates(x), coordinates(y))
	t1 = as.numeric(index(x)) # time in seconds
	t2 = as.numeric(index(y)) # time in seconds
	dt = abs(outer(t1, t2, "-"))
	# separable, product covariance model:
	0.6 * exp(-du/750000) * exp(-dt / (1.5 * 3600 * 24))
}

n = 10
tgrd = seq(min(index(w)), max(index(w)), length=n)
pred = krige0(sqrt(values)~1, w, STF(grd, tgrd), covfn)
layout = list(list("sp.points", pts, first=F, cex=.5),
	list("sp.lines", m, col='grey'))
wind.pr0 = STFDF(grd, tgrd, data.frame(var1.pred = pred))

v = vgmST("separable",
          space = vgm(1, "Exp", 750000), 
          time = vgm(1, "Exp", 1.5 * 3600 * 24),
          sill = 0.6)
wind.ST = krigeST(sqrt(values)~1, w, STF(grd, tgrd), v)

all.equal(wind.pr0, wind.ST)
#}

# stars:
df = data.frame(a = rep(NA, 324*10))
s = STF(grd, tgrd)
newd = addAttrToGeom(s, df)
wind.sta = krigeST(sqrt(values)~1, st_as_stars(w), st_as_stars(newd), v)
# 1
plot(stars::st_as_stars(wind.ST), breaks = "equal", col = sf.colors())
# 2
stplot(wind.ST)
# 3
plot(wind.sta, breaks = "equal", col = sf.colors())
st_as_stars(wind.ST)[[1]][1:3,1:3,1]
(wind.sta)[[1]][1:3,1:3,1]
st_bbox(wind.sta)
bbox(wind.ST)
all.equal(wind.sta, stars::st_as_stars(wind.ST), check.attributes = FALSE)

# 4: roundtrip wind.sta->STFDF->stars
rt = stars::st_as_stars(as(wind.sta, "STFDF"))
plot(rt, breaks = "equal", col = sf.colors())
# 5:
stplot(as(wind.sta, "STFDF"))
st_bbox(rt)

# 6:
stplot(as(st_as_stars(wind.ST), "STFDF"))
