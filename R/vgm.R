# $Id: vgm.q,v 1.9 2006-02-10 19:01:07 edzer Exp $

"vgm" <-
function(psill = 0, model, range = 0, nugget, add.to, anis, kappa = 0.5,
		..., covtable) {
	add.to.df = function(x, y) {
		x = rbind(x, y)
		row.names(x) = 1:nrow(x)
		return(x)
	}
	m = .Call("gstat_get_variogram_models", as.integer(0))
	n = length(m)
	mf = factor(m, levels = m)
	if (missing(model)) {
		ml = .Call("gstat_get_variogram_models", as.integer(1))
		mlf = factor(ml, levels = ml)
		return(data.frame(short = mf, long = mlf))
	}
	table = NULL
	if (model == "Tab" && !missing(covtable)) {
		table = as.matrix(covtable)
		if (NCOL(table) != 2)
			stop("covtable should be a 2-column matrix with distance and cov.")
		range = max(table[,1])
		if (min(table[,1]) != 0.0)
			stop("the first covariance value should be at distance 0.0")
		table = table[,2]
		mf = factor(c(m, "Tab"), levels = c(m, "Tab"))
		if (!missing(add.to) || !missing(nugget))
			stop("cannot add submodels or nugget to covariance Table model")
	} else if (!any(m == model)) 
		stop(paste("variogram model", model, "unknown\n"))
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
	if (!is.null(table))
		attr(ret, "table") = table
	class(ret) = c("variogramModel", "data.frame")
	ret
}
