data(meuse)
data(meuse.grid)
vgm.fit = fit.variogram(variogram(zinc~1, ~x+y, meuse), vgm(1, "Sph", 800, 1))
bl0 = krige(zinc~1, ~x+y, meuse, meuse.grid, model = vgm.fit, block = c(0,0))
bl1 = krige(zinc~1, ~x+y, meuse, meuse.grid, model = vgm.fit, block = c(40,40))
bl2 = krige(zinc~1, ~x+y, meuse, meuse.grid, model = vgm.fit,block = c(100,100))
bl3 = krige(zinc~1, ~x+y, meuse, meuse.grid, model = vgm.fit,block = c(400,400))
bl0$"block=0x0" =     bl0$var1.pred
bl0$"block=40x40" =   bl1$var1.pred
bl0$"block=100x100" = bl2$var1.pred
bl0$"block=400x400" = bl3$var1.pred

library(lattice)

plt1 = levelplot(z ~ x + y | name, map.to.lev(bl0, z = c(5:8)), layout=c(4,1),
        asp=mapasp(bl0),col.regions=bpy.colors(), main = "kriging predictions")
bl0$"block=0x0" =     bl0$var1.var
bl0$"block=40x40" =   bl1$var1.var
bl0$"block=100x100" = bl2$var1.var
bl0$"block=400x400" = bl3$var1.var
plt2 = levelplot(sqrt(z) ~ x + y | name, map.to.lev(bl0, z = c(5:8)), 
	layout=c(4,1), asp=mapasp(bl0),col.regions=bpy.colors(), 
	main = "kriging standard errors")
print(plt1, split = c(1, 1, 1, 2), more = T)
print(plt2, split = c(1, 2, 1, 2), more = F)
# block krige the full area:
bl = krige(zinc~1, ~x+y, meuse, newdata = data.frame(x=0,y=0),
        model = vgm.fit, block = meuse.grid[c("x","y")])
bl
# block kriging standard error:
sqrt(bl$var1.var)
# classical statistical standard error of mean:
sqrt(var(meuse$zinc)/155)
