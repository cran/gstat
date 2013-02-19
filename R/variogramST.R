VgmFillNA = function(x, boundaries) {
  # pads the sample variogram with NA rows where no data are available.
	n = length(boundaries) - 1
	ix = rep(NA, n)
	#ix[which(1:n %in% findInterval(x$dist, boundaries))] = 1:nrow(x)
	ix[ findInterval(x$dist, boundaries)] = 1:nrow(x)
	# x$b = boundaries[-1]
	x[ix,]
}

VgmAverage = function(ret, boundaries) {
	# take out NULL variograms:
	ret = ret[!sapply(ret, is.null)]
	# take care of missing rows...
	ret = lapply(ret, VgmFillNA, 
			boundaries = c(0, 1e-6 * boundaries[2], boundaries[-1]))
	# weighted average:
  # sum np, weighted sum of gamma and dist devided by summed np
	np = apply(do.call(cbind, lapply(ret, function(x) x$np)), 1, sum,
	           na.rm = TRUE)
	gamma = apply(do.call(cbind, lapply(ret, function(x) x$gamma*x$np)), 1, sum,
		na.rm = TRUE)/np
	dist = apply(do.call(cbind, lapply(ret, function(x) x$dist*x$np)), 1, sum,
		na.rm = TRUE)/np
	v = data.frame(np = np, dist = dist, gamma = gamma)
	class(v) = class(ret[[1]])
	attr(v, "boundaries") = attr(ret[[1]], "boundaries")
	v[is.na(v)] = NA
	v
}

StVgmLag = function(formula, data, dt, pseudo, ...) {
	.ValidObs = function(formula, data)
		!is.na(data[[as.character(as.list(formula)[[2]])]])
	d = dim(data)
	ret = vector("list", d[2] - dt)
	if (dt == 0) {
		for (i in 1:d[2]) {
			d0 = data[,i]
      valid = .ValidObs(formula, d0)
      if(sum(valid) <= 1)
        ret[[i]] <- NULL
      else {
        d0 = d0[valid,]
        ret[[i]] = variogram(formula, d0, ...)
      }
		}
	} else {
		for (i in 1:(d[2] - dt)) {
			d1 = data[, i]
      valid1 = .ValidObs(formula, d1)
			d2 = data[, i + dt]
			valid2 = .ValidObs(formula, d2)
      if(sum(valid1)==0 || sum(valid2)==0)
        ret[[i]] <- NULL
      else {
        d1 = d1[valid1,]
        d2 = d2[valid2,]
			  obj = gstat(NULL, paste("D", i, sep=""), formula, d1, 
				  set = list(zero_dist = 3))
			  obj = gstat(obj, paste("D", i+dt, sep=""), formula, d2)
			  ret[[i]] = variogram(obj, cross = "ONLY", pseudo = pseudo, ...)
      }
		}
	}
	VgmAverage(ret, ...)
}

variogramST = function(formula, locations, data, ..., tlags = 0:15, cutoff, 
                       width = cutoff/15, boundaries=seq(0,cutoff,width),
                       progress = TRUE, pseudo = TRUE) {
  if (missing(data))
    data = locations
  if(missing(cutoff))
    cutoff <- spDists(t(data@sp@bbox),!is.projected(data@sp))[1,2]/3
	stopifnot(is(data, "STFDF") || is(data, "STSDF"))
	it = index(data@time)
# 	if (is.regular(
# 				as.zoo(matrix(1:length(it)), order.by = it), 
# 				strict = TRUE)) {
		twidth = diff(it)[1]
		tlags = tlags[tlags <= min(max(tlags), length(unique(it)) - 1)]
# 	} else {
# 		warning("strictly irregular time steps were assumed to be regular")
# 		twidth = mean(diff(it))
# 	}
	ret = vector("list", length(tlags))
	obj = NULL
	t = twidth * tlags
	if (progress)
		pb = txtProgressBar(style = 3, max = length(tlags))
	for (dt in seq(along = tlags)) {
		ret[[dt]] = StVgmLag(formula, data, tlags[dt], pseudo = pseudo, 
                         boundaries = boundaries, ...)
		ret[[dt]]$id = paste("lag", dt - 1, sep="")
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
	b = c(0, b[2]/1e6, b[-1])
	ix = findInterval(v$dist, b)
	b = b[-2]
	spacelags = c(0, b[-length(b)] + diff(b)/2)
	v$spacelag = spacelags[ix]
	if (isTRUE(!is.projected(data)))
		attr(v$spacelag, "units") = "km"
	class(v) = c("StVariogram", "data.frame")
	na.omit(v)
}

plot.StVariogram = function(x, model=NULL, ..., col = bpy.colors(), xlab, ylab, map = TRUE,
		convertMonths = FALSE, wireframe = FALSE, both = FALSE) {
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
  if(!is.null(model)) {
    x[["model"]] <- variogramSurface(model, x[,c("spacelag","timelag")])$model
  }
	x0 = x # needed by wireframe()
	if (!is.null(x$model)) {
		v0 = rbind(x, x)
		v0$what = c(rep("sample", nrow(x)), rep("model", nrow(x)))
		v0$gamma = c(x$gamma, x$model)
		x = v0
	}
	if (wireframe) { 
		if (!is.null(x$model)) {
			if (both)
				wireframe(model+gamma ~ spacelag*timelag, 
					x0, drape = TRUE, col.regions = col, 
					xlab = xlab, ylab = ylab, ...)
			else
				wireframe(model ~ spacelag*timelag, x0, drape = TRUE, 
					col.regions = col, xlab = xlab, ylab = ylab, ...)
		} else
			wireframe(gamma ~ spacelag * timelag, x0, drape = TRUE, col = col,
				xlab = xlab, ylab = ylab, ...)
	} else if (map) {
		if (!is.null(x$model))
			f = gamma ~ spacelag + timelag | what
		else
			f = gamma ~ spacelag + timelag
		levelplot(f, x, xlab = xlab, ylab = ylab, col.regions = col, ...)
	} else { # not map, not wireplot
		if (!is.null(x$model))
			f = gamma ~ dist | what
		else
			f = gamma ~ dist
		x$id = factor(x$id, levels=unique(x$id))
		bp = bpy.colors(length(levels(x$id)))
		ps = list(superpose.line=list(col=bp), superpose.symbol=list(col=bp))
		ylim = c(0, max(x$gamma) * 1.04)
		xlim = c(0, max(x$dist) * 1.04)
		xyplot(f, x, groups = x$id, type='b', ylim = ylim, xlim = xlim,
				auto.key = list(space = "right"), xlab = xlab, 
				par.settings = ps, ...)
	}
}
