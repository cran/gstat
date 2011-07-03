VgmAverage = function(ret) {
	# take care of missing rows...
	gamma = apply(do.call(cbind, lapply(ret, function(x) x$gamma)), 1, mean,
		na.rm = TRUE)
	dist = apply(do.call(cbind, lapply(ret, function(x) x$dist)), 1, mean,
		na.rm = TRUE)
	np = apply(do.call(cbind, lapply(ret, function(x) x$np)), 1, sum,
		na.rm = TRUE)
	v = data.frame(np = np, dist = dist, gamma = gamma)
	class(v) = class(ret[[1]])
	attr(v, "boundaries") = attr(ret[[1]], "boundaries")
	v
}

StVgmLag = function(formula, data, dt, pseudo, ...) {
	d = dim(data)
	ret = vector("list", d[2] - dt)
	if (dt == 0) {
		for (i in 1:d[2])
			ret[[i]] = variogram(formula, data[,i],...)
	} else {
		for (i in 1:(d[2] - dt)) {
			obj = gstat(NULL, paste("D", i, sep=""), formula, data[,i])
			obj = gstat(obj, paste("D", i+dt, sep=""), formula, data[,i+dt])
			ret[[i]] = variogram(obj, cross = "ONLY", pseudo = pseudo, ...)
		}
	}
	VgmAverage(ret)
}

variogramST = function(formula, locations, data, ..., tlags = 0:15,
		progress = TRUE, pseudo = TRUE) {
	if (missing(data))
		data = locations
	stopifnot(is(data, "STFDF") || is(data, "STSDF"))
	it = index(data@time)
	if (is.regular(zoo(1:length(it), it), strict = TRUE)) {
		twidth = diff(it)[1]
		tlags = tlags[tlags <= min(max(tlags), length(unique(it)) - 1)]
	} else {
		warning("strictly irregular time steps were assumed to be regular")
		twidth = mean(diff(it))
	}
	ret = vector("list", length(tlags))
	obj = NULL
	t = twidth * tlags
	if (progress)
		pb = txtProgressBar(style = 3, max = length(tlags))
	for (dt in seq(along = tlags)) {
		ret[[dt]] = StVgmLag(formula, data, tlags[dt], pseudo = pseudo, ...)
		ret[[dt]]$id = paste("lag", dt, sep="")
		if (progress)
			setTxtProgressBar(pb, dt)
	}
	if (progress)
		close(pb)
	# add time lag:
	v = do.call(rbind, ret)
	v$timelag = rep(t, sapply(ret, nrow))
	if (is(t, "yearmon"))
		class(v$timelag) = "yearmon"
	b = attr(ret[[2]], "boundaries")
	ix = findInterval(v$dist, b)
	spacelags = b[1:(length(b)-1)] + diff(b)/2
	v$spacelag = spacelags[ix]
	if (isTRUE(!is.projected(data)))
		attr(v$spacelag, "units") = "km"
	class(v) = c("StVariogram", "data.frame")
	v
}

plot.StVariogram = function(x, ..., col = bpy.colors(), xlab, ylab, map = TRUE,
		convertMonths = FALSE) {
	lst = list(...)
	if (!is.null(lst$col.regions))
		col = lst$col.regions
	if (is(x$timelag, "yearmon")) {
		if (convertMonths) {
			x$timelag = as.numeric(x$timelag) * 12
			attr(x$timelag, "units") = "months"
		} else
			attr(x$timelag, "units") = "years"
	}
	if (missing(xlab)) {
		xlab = "distance"
		u =  attr(x$spacelag, "units")
		if (!is.null(u))
			xlab = paste(xlab, " (", u, ")", sep="")
	}
	if (missing(ylab)) {
		ylab = "time lag"
		u = attr(x$timelag, "units")
		if (!is.null(u))
			ylab = paste(ylab, " (", u, ")", sep="")
	}
	if (!is.null(x$model)) {
		v0 = rbind(x, x)
		v0$what = c(rep("sample", nrow(x)), rep("model", nrow(x)))
		v0$gamma = c(x$gamma, x$model)
		x = v0
	}
	if (map) {
		if (!is.null(x$model))
			f = gamma ~ spacelag + timelag | what
		else
			f = gamma ~ spacelag + timelag
		levelplot(f, x, xlab = xlab, ylab = ylab, col.regions = col, ...)
	} else {
		if (!is.null(x$model))
			f = gamma ~ dist | what
		else
			f = gamma ~ dist
		xyplot(f, x, groups = x$id, type='b', 
				auto.key = list(space = "right"), xlab = xlab, ...)
	}
}

variogram.ST00 = function(formula, locations, data, ...,
		nt = 15, twidth, tcutoff, dX) {
# Sat Jul 16 19:51:37 CEST 2011
	if (missing(data))
		data = locations
	stopifnot(is(data, "STIDF"))
	# da = as(data, "STIDF")
	it = index(data@time)
	if (is.regular(it, strict = TRUE)) {
		if (missing(twidth))
			twidth = diff(it)[1]
		nt = min(nt, length(unique(it))/2)
	} else {
		if (!(missing(twidth) && missing(tcutoff)))
			stop("time is not strictly regular: specify twidth and/or tcutoff")	
		if (!missing(twidth) && !missing(tcutoff)) {
			nt = tcutoff %/% twidth
			if (tcutoff %% twidth > 0)
				nt = nt+1
		}
	}
	if (missing(tcutoff) && !missing(twidth))
		tcutoff = nt * twidth
	if (missing(twidth) && !missing(tcutoff))
		twidth = tcutoff / nt
	# now we have all three; continue with twidth and tcutoff.
	#ret = vector("list", 2*nt+1)
	ret = vector("list", nt)
	obj = NULL
	da = data
	df = as(da, "data.frame")
	dsp = addAttrToGeom(da@sp, df) # now Spatial, with time as attribute
	dsp$RES = residuals(lm(formula, df))
	t = twidth * (0:(nt-1))
	for (dt in 1:nt) {
		id = paste("B", dt-1, sep="")
		df0 = as(dsp, "SpatialPointsDataFrame")
		df0$TIME = df0$time - t[dt] # Back i-1 steps
		obj = gstat(obj, id, RES~TIME, df0)
	}
	if (missing(dX))
		dX = twidth/2
	v = variogram(obj, dX = dX, cross = "ST", progress = TRUE, ...)
	# add time lag:
	tlag = substring(as.character(v$id), first = 5)
	tlag[tlag == ""] = "0"
	v$timelag = t[as.numeric(tlag)+1]
	b = attr(v, "boundaries")
	ix = findInterval(v$dist, b)
	spacelags = b[1:(length(b)-1)] + diff(b)/2
	v$spacelag = spacelags[ix]
	if (isTRUE(!is.projected(data)))
		attr(v$spacelag, "units") = "km"
	class(v) = c("StVariogram", "data.frame")
	v
}
