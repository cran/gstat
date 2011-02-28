# $Id: variogram.default.q,v 1.29 2009-11-02 21:33:17 edzer Exp $

"variogram.default" <-
function (object, locations, X, cutoff, width = cutoff/15.0, alpha = 0, 
    beta = 0, tol.hor = 90/length(alpha), tol.ver = 90/length(beta), 
    cressie = FALSE, dX = numeric(0), boundaries = numeric(0), 
    cloud = FALSE, trend.beta = NULL, debug.level = 1, cross = TRUE, 
	grid, map = FALSE, g = NULL, ..., projected = TRUE, lambda = 1.0) 
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
		cells.dim = length(seq(-cutoff, cutoff, by = width))
		grid.topology = GridTopology(rep(-cutoff, 2), rep(width, 2), rep(cells.dim, 2))
		map = SpatialGrid(grid = grid.topology)
	}
	if (any(is.na(boundaries)))
		stop("no NA values allowed in boundaries")
    .Call("gstat_init", as.integer(debug.level))
	if (!is.null(g) && !is.null(g$set))
		gstat.load.set(g$set)
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
			if (missing(grid) || !is.list(grid)) {
				if (!is.null(g) && gridded(g$data[i]$data))
					grd = unlist(gridparameters(g$data[i]$data))
				else 
					grd = numeric(0)
			} else
				grd = grid[[i]]
            .Call("gstat_new_data", as.double(object[[i]]), 
				as.double(locations[[i]]), as.double(Xloc), 
				as.integer(1), as.double(t.beta), as.integer(-1),
				as.integer(0), as.double(-1), as.integer(1), 
				double(0), grd, as.integer(0), as.integer(projected),
				as.integer(0), as.double(lambda))
			if (!is.null(g) && !is.null(g$model[[id.names[i]]])) 
				load.variogram.model(g$model[[id.names[i]]], c(i - 1, i - 1))
        }
    } else 
		stop("argument object and locations should be lists")

	if (inherits(map, "SpatialGrid"))
		map = as.double(unlist(gridparameters(as(map, "SpatialGrid"))))

    pos = 0
    ids = NULL
	is.direct = NULL
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
                    	as.numeric(dX), as.numeric(boundaries), map)
				  boundaries = numeric(0)
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
					  is.direct = c(is.direct, id1 == id2)
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
							ret.call[[4]], ret.call[[3]])
					  names = c("dx", "dy", id, paste("np", id, sep="."))
					} else { 
					  ret = data.frame(ret, ret.call[[4]], ret.call[[3]])
					  names = c(names, id, paste("np", id, sep="."))
					}
					names(ret) = names
				  }
                }
            }
        }
    }
    .Call("gstat_exit", NULL)
	if (is.logical(map) && map == FALSE) {
    	if (!is.null(ids)) {
			ret$id = factor(ids, levels = unique(ids))
			attr(ret, "direct") = data.frame(id = unique(ids), is.direct = is.direct)
    		if (cloud) {
        		class(ret) = c("variogramCloud", "data.frame")
				attr(ret, ".BigInt") = 2^(4 * .Machine$sizeof.long)
    		} else 
				class(ret) = c("gstatVariogram", "data.frame")
		}
	} else {
		coordinates(ret) = c("dx", "dy")
		gridded(ret) = TRUE
		ret = list(map = ret)
		class(ret) = c("variogramMap", "list")
	}
    ret
}
