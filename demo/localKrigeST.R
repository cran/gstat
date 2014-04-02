## FNN local prediction
########################
library(sp)
library(spacetime)
library(gstat)

# create n space-time points over [0,1] x [0,1] x [Now, Now+some days]
t0 = Sys.time() # now
n = 1000
set.seed(13131) # fix outcomes
x = runif(n)
y = runif(n)
t = t0 + 1e6 * runif(n)
z = rnorm(n)
stidf = STIDF(SpatialPoints(cbind(x,y)),sort(t),data.frame(z=z))

# create a regular 20 x 20 x 10 grid of prediction locations:
grd = as(SpatialGrid(GridTopology(c(0.025,0.025), c(.05, .05), c(20,20))), "SpatialPixels")
tgrd = seq(min(t), max(t), length.out = 10)
stf = STF(grd, tgrd)

# define a variogram model
sumMetricModel <- vgmST("sumMetric",
                        space=vgm(0, "Exp", 1),
                        time =vgm(0, "Exp",  1),
                        joint=vgm(0.8, "Exp", 1, 0.2),
                        stAni=1/1e6)
attr(sumMetricModel, "temporal unit") <- "secs"

locKrig <- krigeST(z~1, stidf, stf, sumMetricModel, nmax=50)
stplot(locKrig, col.regions=bpy.colors(), scales=list(draw=T))