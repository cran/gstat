data(meuse)
data(meuse.grid)

# ordinary kriging
v <- variogram(log(zinc)~1,~x+y, meuse)
m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
plot(v, model = m)
lzn.kr <- krige(formula = log(zinc)~1, locations = ~x+y, model = m, 
	data = meuse, newdata = meuse.grid)
pl1 <- levelplot(var1.pred ~ x + y, lzn.kr, aspect = mapasp(lzn.kr),
	main = "ordinary kriging prediction of log-zinc")
pl2 <- levelplot(sqrt(var1.var) ~ x + y, lzn.kr, aspect = mapasp(lzn.kr),
	main = "ordinary kriging prediction error")

# obtain ndist values in meuse.grid at data locations of meuse:
# (inverse distance is fine, as the neighbourhood size is 1)
x <- krige(ndist ~ 1, ~ x + y, meuse.grid, meuse, model = NULL, nmax = 1)
meuse$ndist = x$var1.pred

# universal kriging
v <- variogram(log(zinc)~sqrt(ndist),~x+y, meuse)
m <- fit.variogram(v, vgm(1, "Exp", 300, 1))
plot(v, model = m)
lzn.kr <- krige(formula = log(zinc)~sqrt(ndist), locations = ~x+y, model = m, 
	data = meuse, newdata = meuse.grid)
pl3 <- levelplot(var1.pred ~ x + y, lzn.kr, aspect = mapasp(lzn.kr),
	main = "universal kriging prediction of log-zinc")
pl4 <- levelplot(sqrt(var1.var) ~ x + y, lzn.kr, aspect = mapasp(lzn.kr),
	main = "universal kriging prediction error")
print(pl1, split = c(1,1,2,2), more = T)
print(pl2, split = c(1,2,2,2), more = T)
print(pl3, split = c(2,1,2,2), more = T)
print(pl4, split = c(2,2,2,2))
