"vgm" <-
function(psill = 0, model, range = 0, nugget, add.to, anis, kappa = 0.5) {
	add.to.df = function(x, y) {
		x = rbind(y, x)
		row.names(x) = 1:nrow(x)
		return(x)
	}
	n = .Call("gstat_get_n_variogram_models", 0, PACKAGE = "gstat")[[1]];
	m = .C("Cgstat_get_variogram_models", rep("",n), PACKAGE = "gstat")[[1]]
	mf = factor(m, levels = m)
	if (missing(model))
		return(mf)
	if (!any(m == model)) stop(paste("variogram model", model, "unknown\n"))
	if (missing(anis))
		anis = c(0,0,0,1,1)
	if (length(anis) == 2)
		anis = c(anis[1], 0, 0, anis[2], 1)
	else if (length(anis) != 5)
		stop("anis vector should have length 2 (2D) or 5 (3D)")
	if (model != "Nug") {
		if (model != "Lin" && model != "Err" && model != "Int")
			if (range <= 0.0) stop("range should be positive")
		else if(range < 0.0) stop("range should be non-negative")
	} else {
		if (range != 0.0) stop("Nugget should have zero range")
		if (anis[4] != 1.0 || anis[5] != 1.0)
			stop("Nugget anisotropy is nonsense")
	}
	if (!missing(nugget)) {
		ret = data.frame(model=mf[mf==model], psill=psill, range=range,
			kappa = kappa, ang1=anis[1], ang2=anis[2], ang3=anis[3], 
			anis1=anis[4], anis2=anis[5])
		n.vgm = data.frame(model=mf[mf=="Nug"], psill=nugget, range=0,
			kappa = 0.0, ang1=0.0, ang2=0.0, ang3=0.0, anis1=1.0, anis2=1.0)
		ret = add.to.df(n.vgm, ret)
	} else
		ret = data.frame(model=mf[mf==model], psill=psill, range=range,
			kappa = kappa, ang1=anis[1], ang2=anis[2], ang3=anis[3], 
			anis1=anis[4], anis2=anis[5])
	if (!missing(add.to))
		ret = add.to.df(data.frame(add.to), ret)
	class(ret) = c("variogramModel", "data.frame")
	ret
}
