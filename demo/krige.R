# $Id: krige.R,v 1.4 2006-02-10 19:05:02 edzer Exp $
data(meuse)
data(meuse.grid)

# ordinary kriging
v <- variogram(log(zinc)~1,~x+y, meuse)
m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
plot(v, model = m)
lzn.kr <- krige(formula = log(zinc)~1, locations = ~x+y, model = m, 
	data = meuse, newdata = meuse.grid)

library(lattice)

pl1 <- levelplot(var1.pred ~ x + y, lzn.kr, aspect = "iso",
	main = "ordinary kriging prediction of log-zinc")
pl2 <- levelplot(sqrt(var1.var) ~ x + y, lzn.kr, aspect = "iso",
	main = "ordinary kriging prediction error")

# universal kriging
v <- variogram(log(zinc)~sqrt(dist),~x+y, meuse)
m <- fit.variogram(v, vgm(1, "Exp", 300, 1))
plot(v, model = m)
lzn.kr <- krige(formula = log(zinc)~sqrt(dist), locations = ~x+y, model = m, 
	data = meuse, newdata = meuse.grid)
pl3 <- levelplot(var1.pred ~ x + y, lzn.kr, aspect = "iso",
	main = "universal kriging prediction of log-zinc")
pl4 <- levelplot(sqrt(var1.var) ~ x + y, lzn.kr, aspect = "iso",
	main = "universal kriging prediction error")
print(pl1, split = c(1,1,2,2), more = T)
print(pl2, split = c(1,2,2,2), more = T)
print(pl3, split = c(2,1,2,2), more = T)
print(pl4, split = c(2,2,2,2))
