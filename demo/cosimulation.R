data(meuse)
data(meuse.grid)

# cosimulation the four heavy metal variables
meuse.g <- gstat(id="zn", formula=zinc~1,loc=~x+y, data=meuse, nmax = 10,
	set = list(zero = 1e-10))
meuse.g <- gstat(meuse.g, "cu", copper~1, ~x+y, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cd", cadmium~1,~x+y, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "pb", lead~1,   ~x+y, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.g, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.g)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid, nsim = 2)
pl1 <- levelplot(z ~ x + y | name, map.to.lev(z, z=c(3:4)), aspect = mapasp(z),
	main = "zinc simulations")
pl2 <- levelplot(z ~ x + y | name, map.to.lev(z, z=c(5:6)), aspect = mapasp(z),
	main = "copper simulations")
pl3 <- levelplot(z ~ x + y | name, map.to.lev(z, z=c(7:8)), aspect = mapasp(z),
	main = "cadmium simulations")
pl4 <- levelplot(z ~ x + y | name, map.to.lev(z, z=c(9:10)), aspect = mapasp(z),
	main = "lead simulations")
print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))

# indicator cosimulation for the 9 percentiles of zinc:
data(meuse)
data(meuse.grid)
q <- quantile(meuse$zinc, seq(.1,.9,.1))
meuse.i <- gstat(id = "zn1", formula = I(zinc < q[1])~1, loc = ~x+y, 
	data = meuse, nmax = 7, beta = .1, set = list(order = 4, zero = 1e-5))
meuse.i <- gstat(meuse.i, "zn2", I(zinc < q[2])~1, ~x+y, meuse, 
	nmax = 7, beta=.2)
meuse.i <- gstat(meuse.i, "zn3", I(zinc < q[3])~1, ~x+y, meuse, 
	nmax = 7, beta=.3)
meuse.i <- gstat(meuse.i, "zn4", I(zinc < q[4])~1, ~x+y, meuse, 
	nmax = 7, beta=.4)
meuse.i <- gstat(meuse.i, "zn5", I(zinc < q[5])~1, ~x+y, meuse, 
	nmax = 7, beta=.5)
meuse.i <- gstat(meuse.i, "zn6", I(zinc < q[6])~1, ~x+y, meuse, 
	nmax = 7, beta=.6)
meuse.i <- gstat(meuse.i, "zn7", I(zinc < q[7])~1, ~x+y, meuse, 
	nmax = 7, beta=.7)
meuse.i <- gstat(meuse.i, "zn8", I(zinc < q[8])~1, ~x+y, meuse, 
	nmax = 7, beta=.8)
meuse.i <- gstat(meuse.i, "zn9", I(zinc < q[9])~1, ~x+y, meuse, 
	nmax = 7, beta=.9)

meuse.i <- gstat(meuse.i, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.i, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.i)
plot(x, model = meuse.fit)

z <- predict(meuse.fit, newdata = meuse.grid, nsim = -2)
# nsim = -2 instead of 2 forces indicator simulation
levelplot(z~x+y|name, map.to.lev(z, z=c(3:20)), aspect = mapasp(z))
rm(z, meuse.fit, x, meuse.i, q)
