## ex01.cmd, ex02.cmd:
##
## Two variables with (initial estimates of) variograms,
## calcute sample variogram and plot fitted model
##
par(ask = TRUE)
data(meuse)
x <- variogram(zinc ~ 1, ~ x + y, meuse)
v <- vgm(140000, "Sph", 800, nug = 10000)
plot(x, model = v)
plot(x, model = fit.variogram(x, model = v))
x <- variogram(log(zinc) ~ 1, ~ x + y, meuse)
v <- vgm(.5, "Sph", 800, nug = .1)
plot(x, model = v)
plot(x, model = fit.variogram(x, model = v))
##
## ex03.cmd:
## Inverse distance interpolation on a mask map
##
data(meuse.grid)
x <- krige(zinc ~ 1, ~ x + y, meuse, meuse.grid, model = NULL)
levelplot(var1.pred ~ x + y, x, aspect = mapasp(x))
##
## ex04.cmd 
## Local ordinary block kriging at non-gridded locations
##
## skipped min and radius data arguments, which affects only the
## last location
new.locs <- data.frame(x = c(181170, 180310, 180205, 178673, 178770, 178270),
	y = c(333250, 332189, 331707, 330066, 330675, 331075))
krige(zinc ~ 1, ~ x + y, meuse, newdata = new.locs, 
		model = vgm(1.34e5, "Sph", 800, nug = 2.42e4), 
		block = c(40,40), nmax = 40)

##
## ex05.cmd 
##
## Local simple point kriging on a mask map
##
v <- vgm(0.581, "Sph", 900, nug = 0.0554)
x <- krige(log(zinc) ~ 1, ~ x + y, meuse, meuse.grid, model = v, 
	nmax = 40, beta = c(5.9))
levelplot(var1.pred ~ x + y, x, aspect = mapasp(x), 
	main = "log(zinc) simple kriging prediction")
levelplot(sqrt(var1.var) ~ x + y, x, aspect = mapasp(x), 
	main = "log(zinc) simple kriging standard errors")
##
## ex06.cmd 
##
## Unconditional Gaussian simulation on a mask
## (local neigbourhoods, simple kriging)
##
x <- krige(log(zinc) ~ 1, ~ x + y, data = NULL, newdata = meuse.grid, 
	model = v, nmax = 20, beta = c(5.9), nsim = 5, dummy = TRUE)
levelplot(z ~ x + y | name, map.to.lev(x, z=c(3:7)), aspect = mapasp(x),
	main = "five unconditional realisations of a correlated Gaussian field")
##
## ex07.cmd 
##
## Gaussian simulation, conditional upon data
## (local neighbourhoods, simple and ordinary kriging)
##
x <- krige(log(zinc) ~ 1, ~ x + y, meuse, meuse.grid, 
	model = v, nmax = 20, beta = c(5.9), nsim = 5)
levelplot(z ~ x + y | name, map.to.lev(x, z=c(3:7)), aspect = mapasp(x),
	main = "five conditional realisations of a correlated Gaussian field")

##
## ex08.cmd 
##
## Change of support: local ordinary block kriging on a mask
##
x <- krige(log(zinc) ~ 1, ~ x + y, meuse, meuse.grid, 
	model = v, nmax = 40, block = c(40,40))
levelplot(var1.pred ~ x + y, x, aspect = mapasp(x),
	main = "ordinary block kriging predictions")
levelplot(sqrt(var1.var) ~ x + y, x, aspect = mapasp(x),
	main = "ordinary block kriging prediction standard errors")

##
## ex09.cmd 
##
## Obtain map values at data() locations
## (Point-map overlay)
##
# we trick here by using inv.weighted distance interpolation, using the
# single nearest observation. It will not fail on points outside the grid.
# Note that we reversed meuse.grid and meuse to get these results.
# also, I prefer not to use the sqrt(ndist) variable, as S can do sqrt()
x <- krige(ndist ~ 1, ~ x + y, meuse.grid, meuse, model = NULL, nmax = 1)
plot(meuse$dist ,x$var1.pred) # interesting! (note the difference in units)
# use the new ones:
meuse$ndist = x$var1.pred
x <- krige(part.a ~ 1, ~ x + y, meuse.grid, meuse, model = NULL, nmax = 1)
meuse$part.a = x$var1.pred
x <- krige(part.b ~ 1, ~ x + y, meuse.grid, meuse, model = NULL, nmax = 1)
meuse$part.b = x$var1.pred

##
## ex10.cmd 
##
## Multiple kriging: prediction on more than one variable
## (ordinary kriging of two variables)
## (note that zinc_map.eas wass obtained through ex09.gst)
##
x <- variogram(ndist~1,~x+y,meuse)
v.ndist <- fit.variogram(x, model = vgm(1,"Gau",100))
plot(x, model = v.ndist)
g <- gstat(id = "ln.zinc", form = log(zinc) ~ 1, loc = ~ x + y, 
	data = meuse, nmax = 40, model = v)
g <- gstat(g, id = "ndist", form = ndist ~ 1, loc = ~ x + y, 
	data = meuse, nmax = 40, model = vgm(.01, "Nug", 0, add.to = v.ndist))
# the added nugget variance is necessary to avoid near-singular covariances
x <- predict(g, meuse.grid)
levelplot(ln.zinc.pred ~ x + y, x, aspect = mapasp(x),
	main = "log(zinc) ordinary kriging predictions")
levelplot(sqrt(ln.zinc.var) ~ x + y, x, aspect = mapasp(x),
	main = "log(zinc) ordinary kriging prediction standard errors")
levelplot(ndist.pred ~ x + y, x, aspect = mapasp(x),
	main = "ndist ordinary kriging predictions")
levelplot(sqrt(ndist.var) ~ x + y, x, aspect = mapasp(x),
	main = "ndist ordinary kriging prediction standard errors")

##
## ex11.cmd 
##
## Multivariable kriging: ordinary local cokriging of two variables
## For examples of fitting an LMC: see demo(cokriging)
##
g <- gstat(id = "ln.zinc", form = log(zinc) ~ 1, loc = ~ x + y, 
	data = meuse, nmax = 40, model = vgm(0.581, "Sph", 900, 0.0554))
g <- gstat(g, id = "sq.ndist", form = sqrt(ndist) ~ 1, loc = ~ x + y, 
	data = meuse, nmax = 40, model = vgm(0.0631, "Sph", 900, 0.0001))
g <- gstat(g, id = c("ln.zinc", "sq.ndist"), 
	model = vgm(-0.156, "Sph", 900, 1e-5))
# small nugget necessary to let gstat recognize LMC
x <- predict(g, meuse.grid)
levelplot(ln.zinc.pred ~ x + y, x, aspect = mapasp(x),
	main = "log(zinc) ordinary cokriging predictions")
levelplot(sqrt(ln.zinc.var) ~ x + y, x, aspect = mapasp(x),
	main = "log(zinc) ordinary cokriging prediction standard errors")
levelplot(sq.ndist.pred ~ x + y, x, aspect = mapasp(x),
	main = "ndist ordinary cokriging predictions")
levelplot(sqrt(sq.ndist.var) ~ x + y, x, aspect = mapasp(x),
	main = "ndist ordinary cokriging prediction standard errors")

##
## ex12.cmd 
##
## Stratified ordinary kriging (within-category ordinary kriging)
##
## (stratified mode not implemented)

##
## ex13.cmd 
##
## Local universal kriging, using one continuous variable
###
## the variogram should be that of the residual:
x <- krige(log(zinc) ~ sqrt(ndist), ~ x + y, meuse, meuse.grid, 
	model = vgm(.149, "Sph", 700, .0674), nmax = 40)
levelplot(var1.pred ~ x + y, x, aspect = mapasp(x),
	main = "universal kriging predictions")
levelplot(sqrt(var1.var) ~ x + y, x, aspect = mapasp(x),
	main = "universal kriging prediction standard errors")

##
## ex14.cmd 
##
## Universal kriging, using one continuous and
## two binary variables.
##
x <- krige(log(zinc) ~ -1 + sqrt(ndist)+ part.a + part.b, 
	~ x + y, meuse, meuse.grid, 
	model = vgm(.149, "Sph", 700, .0674))
levelplot(part.a ~ x + y, meuse.grid, aspect = mapasp(meuse.grid),
	main = "the areas defining part.a (1) and part.b (0)")
levelplot(var1.pred ~ x + y, x, aspect = mapasp(x),
	main = "universal kriging predictions")
levelplot(sqrt(var1.var) ~ x + y, x, aspect = mapasp(x),
	main = "universal kriging prediction standard errors")

##
## ex14a.cmd 
## 
## stratified universal kriging: 
## (again: not implemented)
##

## ex15.cmd 
##
## Local linear model, using one continuous variable
##
x <- krige(log(zinc) ~ sqrt(ndist), ~ x + y, meuse, meuse.grid, 
	model = NULL, nmax = 40)
levelplot(var1.pred ~ x + y, x, aspect = mapasp(x),
	main = "IID local linear model kriging predictions")
levelplot(sqrt(var1.var) ~ x + y, x, aspect = mapasp(x),
	main = "IID local linear model prediction standard errors")

##
## ex16.cmd 
##
## Multivariable indicator cosimulation 
## ==>> see demo(cosimulation) for an extended example how to do this
##

##
## ex17.cmd 
##
## global coordinate polynomial trend surfaces
## trend orders 0-3.
## (you'd better use lm() for this)
##
