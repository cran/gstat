# simple demo of 3D interpolation of 50 random points in the unit cube
n <- 50
data3D <- data.frame(x = runif(n), y = runif(n), z = runif(n), v = rnorm(n))
range1D <- seq(from=0,to=1,length=20)
grid3D <- expand.grid(range1D, range1D, range1D)
names(grid3D) <- c("x", "y", "z")
x <- krige(formula = v~1, locations = ~x+y+z, 
	data = data3D, newdata = grid3D, model = vgm(1, "Exp", .2))
levelplot(var1.pred ~ x + y | z, x)
rm(n, data3D, range1D, grid3D, x)
