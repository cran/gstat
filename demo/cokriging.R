data(meuse)
data(meuse.grid)

# cokriging of the four heavy metal variables
meuse.g <- gstat(id="zn", formula=log(zinc)~1,loc=~x+y, data=meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cu", log(copper)~1, ~x+y, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "cd", log(cadmium)~1,~x+y, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, "pb", log(lead)~1,   ~x+y, meuse, nmax = 10)
meuse.g <- gstat(meuse.g, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.g, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.g)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid)
pl1 <- levelplot(zn.pred~x+y, z, aspect = mapasp(z), 
	main="log-zinc predictions")
pl2 <- levelplot(cu.pred~x+y, z, aspect = mapasp(z), 
	main="log-copper predictions")
pl3 <- levelplot(cd.pred~x+y, z, aspect = mapasp(z), 
	main="log-cadmium predictions")
pl4 <- levelplot(pb.pred~x+y, z, aspect = mapasp(z), 
	main="log-lead predictions")
print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))
pl1 <- levelplot(sqrt(zn.var)~x+y, z, aspect = mapasp(z), 
	main="log-zinc std.err.")
pl2 <- levelplot(sqrt(cu.var)~x+y, z, aspect = mapasp(z), 
	main="log-copper std.err.")
pl3 <- levelplot(sqrt(cd.var)~x+y, z, aspect = mapasp(z), 
	main="log-cadmium std.err.")
pl4 <- levelplot(sqrt(pb.var)~x+y, z, aspect = mapasp(z), 
	main="log-lead st.err.")
print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))

rm(meuse.g, x, meuse.fit, z)

# indicator cokriging for the 9 percentiles of zinc:
data(meuse)
data(meuse.grid)
q <- quantile(meuse$zinc, seq(.1,.9,.1))
meuse.i <- gstat(id = "zn1", formula = I(zinc < q[1])~1, loc = ~x+y, 
	data = meuse, nmax = 7, beta = .1, set = list(order = 4, zero = 1e-5))
meuse.i <- gstat(meuse.i, "zn2", I(zinc < q[2])~1, ~x+y, meuse, nmax = 7, beta=.2)
meuse.i <- gstat(meuse.i, "zn3", I(zinc < q[3])~1, ~x+y, meuse, nmax = 7, beta=.3)
meuse.i <- gstat(meuse.i, "zn4", I(zinc < q[4])~1, ~x+y, meuse, nmax = 7, beta=.4)
meuse.i <- gstat(meuse.i, "zn5", I(zinc < q[5])~1, ~x+y, meuse, nmax = 7, beta=.5)
meuse.i <- gstat(meuse.i, "zn6", I(zinc < q[6])~1, ~x+y, meuse, nmax = 7, beta=.6)
meuse.i <- gstat(meuse.i, "zn7", I(zinc < q[7])~1, ~x+y, meuse, nmax = 7, beta=.7)
meuse.i <- gstat(meuse.i, "zn8", I(zinc < q[8])~1, ~x+y, meuse, nmax = 7, beta=.8)
meuse.i <- gstat(meuse.i, "zn9", I(zinc < q[9])~1, ~x+y, meuse, nmax = 7, beta=.9)
meuse.i <- gstat(meuse.i, model=vgm(1, "Sph", 900, 1), fill.all=T)
x <- variogram(meuse.i, cutoff=1000)
meuse.fit = fit.lmc(x, meuse.i)
plot(x, model = meuse.fit)
z <- predict(meuse.fit, newdata = meuse.grid)
levelplot(z ~ x + y | name, map.to.lev(z, z=c(3,5,7,9,11,13,15,17,19),
	 ns = paste("est.Pr(Zn < ", q, ")", sep = "")), aspect = mapasp(z))
rm(z, meuse.fit, x, meuse.i, q)
