"variogram.default" <-
function (object, locations, X, cutoff, width = cutoff/15, alpha = 0, 
    beta = 0, tol.hor = 90/length(alpha), tol.ver = 90/length(beta), 
    cressie = FALSE, dX = numeric(0), boundaries = numeric(0), 
    cloud = FALSE, trend.beta = NULL, debug.level = 1, cross = TRUE, 
	grid, map = FALSE, ...) 
{
    id1 = id2 = 0
    ret = NULL
    if (missing(cutoff)) {
		if (is.logical(map) && map == TRUE)
			stop("for variogram maps, supply at least a cutoff")
        cutoff = numeric(0)
	}
    else if (width <= 0) 
        stop("argument width should be positive")
    if (missing(width)) 
        width = numeric(0)
    if (cloud == TRUE) 
        width = 0
	if (is.logical(map) && map == TRUE) {
		# set up map:
		dx = seq(-cutoff, cutoff, by = width)
		map = data.frame(x = dx, y = dx)
		map = data.frame(lapply(map, as.double))
		require(sp)
		coordinates(map) = c("x", "y")
		gridded(map) = TRUE
	}
    .Call("gstat_init", as.integer(debug.level), PACKAGE = "gstat")
    id.names = NULL
    if (is.list(object) && is.list(locations)) {
        nvars = length(object)
        if (!is.null(names(object))) 
            id.names = names(object)
        for (i in 1:nvars) {
            if (missing(X)) 
                Xloc = rep(1, length(object[[i]]))
            else Xloc = X[[i]]
            t.beta = numeric(0)
            if (!is.null(trend.beta) && length(trend.beta) > 0) 
                t.beta = trend.beta[[i]]
            else t.beta = numeric(0)
			if (missing(grid) || !is.list(grid))
				grd = numeric(0)
			else
				grd = grid[[i]]
            .Call("gstat_new_data", as.double(object[[i]]), 
				as.double(locations[[i]]), as.double(Xloc), 
				as.integer(1), as.double(t.beta), as.integer(-1),
				as.integer(0), as.double(-1), as.integer(1), 
				double(0), grd
				, PACKAGE = "gstat"
				)
        }
    } else 
		stop("argument object and locations should be lists")

	if (is(map, "SpatialDataFrameGrid"))
		map = as.double(unlist(gridparameters(map)))

    pos = 0
    ids = NULL
	if (cross)
		id.range = nvars:1
	else
		id.range = 1:nvars
    for (id1 in id.range) {
        for (id2 in ifelse(cross, 1, id1):id1) {
            if (is.null(id.names)) 
                id = ifelse(id1 == id2, paste(id1), cross.name(id2, id1))
            else id = ifelse(id1 == id2, paste(id.names[id1]), 
                cross.name(id.names[id2], id.names[id1]))
            for (a in alpha) {
                for (b in beta) {
                  direction = as.numeric(c(a, b, tol.hor, tol.ver))
                  ret.call = .Call("gstat_variogram", 
				  		as.integer(c(id1 - 1, id2 - 1)), 
						as.numeric(cutoff), as.numeric(width), 
                    	as.numeric(direction), as.integer(cressie), 
                    	as.numeric(dX), as.numeric(boundaries), map
						, PACKAGE = "gstat"
						)
                  if (is.logical(map) && map == FALSE) {
                    np = ret.call[[1]]
                    sel = np > 0
                    n.dir = length(sel[sel])
                    if (n.dir > 0) {
                      dist = ret.call[[2]]
                      gamma = ret.call[[3]]
                      dir.a = rep(a, n.dir)
                      dir.b = rep(b, n.dir)
                      ids = c(ids, rep(id, n.dir))
                      df = data.frame(np = np[sel], dist = dist[sel], 
                        gamma = gamma[sel], dir.hor = dir.a, dir.ver = dir.b)
                      if (pos > 0) 
                        ret[(pos + 1):(pos + n.dir), ] = df
                      else ret = df
                      pos = pos + n.dir
                    }
				  } else {
				  	if (is.null(ret)) {
					  ret = data.frame(ret.call[[1]], ret.call[[2]],
							ret.call[[3]], ret.call[[4]])
					  names = c("dx", "dy", paste("np", id, sep="."), id)
					} else { 
					  ret = data.frame(ret, ret.call[[3]], ret.call[[4]])
					  names = c(names, paste("np", id, sep="."), id)
					}
					names(ret) = names
				  }
                }
            }
        }
    }
    .Call("gstat_exit", NULL, PACKAGE = "gstat")
	if (is.logical(map) && map == FALSE) {
    	ret$id = factor(ids, levels = unique(ids))
    	if (cloud) 
        	class(ret) = c("variogramCloud", "data.frame")
    	else 
			class(ret) = c("gstatVariogram", "data.frame")
	} else {
		require(sp)
		coordinates(ret) = c("dx", "dy")
		gridded(ret) = TRUE
		ret = list(map = ret)
		class(ret) = c("variogramMap", "list")
	}
    ret
}
