if (require(sp)) {
set.seed(13331)

x = runif(10)
y = runif(10)
z = rnorm(10)
d = data.frame(x=x,y=y,z=z)
bl = c(1,1)

library(gstat)
nd = data.frame(x=.5, y=.5)

# single block:
library(sp)
coordinates(d) = ~x+y
nd = SpatialPoints(nd)
xx = krige(z~1, d, nd, model=vgm(1, "Exp", 1), block = bl,
	set = list(nb=4))

ring = cbind(c(0,bl[1],bl[1],0,0), c(0,0,bl[2],bl[2],0))
r1 = SpatialRings(list(Srings(list(Sring(ring)), ID = "xx")))
a = data.frame(a = 1, b = 2)
rownames(a) = "xx"
r1df = SpatialRingsDataFrame(r1, a)

g = gstat(formula=z~1, data=d, model=vgm(1, "Exp", 1))
args = list(type = "regular", n=16, offset=c(0.5,0.5))
yy = predict(g, r1df, block = bl, sps.args = args)

print(all.equal(as.data.frame(xx), as.data.frame(yy)))

## multiple blocks of equal size:
nd = data.frame(x= 0:4 + .5, y=rep(.5,5))
nd = SpatialPoints(nd)
xx = krige(z~1, d, nd, model=vgm(1, "Exp", 1), block = bl,
	set = list(nb=4))
ring0 = cbind(c(0,bl[1],bl[1],0,0), c(0,0,bl[2],bl[2],0))
ring1 = cbind(c(1+0,1+bl[1],1+bl[1],1+0,1+0), c(0,0,bl[2],bl[2],0))
ring2 = cbind(c(2+0,2+bl[1],2+bl[1],2+0,2+0), c(0,0,bl[2],bl[2],0))
ring3 = cbind(c(3+0,3+bl[1],3+bl[1],3+0,3+0), c(0,0,bl[2],bl[2],0))
ring4 = cbind(c(4+0,4+bl[1],4+bl[1],4+0,4+0), c(0,0,bl[2],bl[2],0))
r1 = SpatialRings(list(
	Srings(list(Sring(ring0)), ID = "x0"),
	Srings(list(Sring(ring1)), ID = "x1"),
	Srings(list(Sring(ring2)), ID = "x2"),
	Srings(list(Sring(ring3)), ID = "x3"),
	Srings(list(Sring(ring4)), ID = "x4"),
	))
df = data.frame(a=rep(1,5), b= rep(1,5))
rownames(df) = c("x0", "x1", "x2", "x3", "x4")
r1df = SpatialRingsDataFrame(r1, df)

yy = predict(g, r1, block = bl, sps.args = args)
all.equal(as.data.frame(xx), as.data.frame(yy))

## multiple blocks of equal size:
args = list(type = "regular", cellsize=.25, offset=c(0.5,0.5))
yy = predict(g, r1, block = bl, sps.args = args)
print(all.equal(as.data.frame(xx), as.data.frame(yy)))

## multiple blocks of varying size:
nd = data.frame(x=c(0.5, 2, 4.5), y=c(0.5, 1, 1.5))
nd = SpatialPoints(nd)
bl = c(1,1)
ring0 = cbind(c(0,bl[1],bl[1],0,0), c(0,0,bl[2],bl[2],0))
xx1 = krige(z~1, d, nd[1], model=vgm(1, "Exp", 1), block = bl,
	set = list(nb=4))
bl = c(2,2)
ring1 = cbind(c(1+0,1+bl[1],1+bl[1],1+0,1+0), c(0,0,bl[2],bl[2],0))
xx2 = krige(z~1, d, nd[2], model=vgm(1, "Exp", 1), block = bl,
	set = list(nb=4))
bl = c(3,3)
ring2 = cbind(c(3+0,3+bl[1],3+bl[1],3+0,3+0), c(0,0,bl[2],bl[2],0))
xx3 = krige(z~1, d, nd[3], model=vgm(1, "Exp", 1), block = bl,
	set = list(nb=4))
r1 = SpatialRings(list(
	Srings(list(Sring(ring0)), ID = "x0"),
	Srings(list(Sring(ring1)), ID = "x1"),
	Srings(list(Sring(ring2)), ID = "x2"),
	))
df = data.frame(a = rep(1,3), b = rep(1,3))
rownames(df) = c("x0", "x1", "x2")
r1df = SpatialRingsDataFrame(r1, df)

args = list(type = "regular", n=16, offset=c(0.5,0.5))
yy = predict(g, r1df, block = bl, sps.args = args)
xx = rbind(as.data.frame(xx1), as.data.frame(xx2), as.data.frame(xx3))
row.names(xx) = 1:3
print(all.equal(as.data.frame(xx), as.data.frame(yy)))
} else
	print("this test only runs if package sp is installed")
