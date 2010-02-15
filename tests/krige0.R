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
x0 <- krige0(zinc~1, meuse, meuse.grid, v, beta = 500, computeVar = TRUE)
x <- krige(zinc~1, meuse, meuse.grid, v, beta = 500)
summary(x$var1.pred - x0$pred)
summary(x$var1.var - x0$var)
# test ok:
x0 <- krige0(zinc~1, meuse, meuse.grid, v, computeVar = TRUE)
x <- krige(zinc~1, meuse, meuse.grid, v)
summary(x$var1.pred - x0$pred)
summary(x$var1.var - x0$var)
# test uk:
x0 <- krige0(zinc~sqrt(dist), meuse, meuse.grid, v, computeVar = TRUE)
x <- krige(zinc~sqrt(dist), meuse, meuse.grid, v)
summary(x$var1.pred - x0$pred)
summary(x$var1.var - x0$var)
