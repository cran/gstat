# test -- load data:
library(gstat)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# test -- idw
meuse.grid$idw <- idw0(zinc~1, meuse, meuse.grid)
x <- idw(zinc~1, meuse, meuse.grid)
summary(x$var1.pred - meuse.grid$idw)
spplot(meuse.grid["idw"],col.regions=bpy.colors())
v = vgm(1, "Exp", 500)
# test sk:
meuse.grid$sk <- krige0(zinc~1, meuse, meuse.grid, v, beta = 500)
x <- krige(zinc~1, meuse, meuse.grid, v, beta = 500)
summary(x$var1.pred - meuse.grid$sk)
spplot(meuse.grid["sk"],col.regions=bpy.colors())
# test ok:
meuse.grid$ok <- krige0(zinc~1, meuse, meuse.grid, v)
x <- krige(zinc~1, meuse, meuse.grid, v)
summary(x$var1.pred - meuse.grid$ok)
# test uk:
meuse.grid$uk <- krige0(zinc~sqrt(dist), meuse, meuse.grid, v)
x <- krige(zinc~sqrt(dist), meuse, meuse.grid, v)
summary(x$var1.pred - meuse.grid$uk)
