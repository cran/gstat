"variogram.default" <-
function (object, locations, X, cutoff, width = cutoff/15, alpha = 0, 
    beta = 0, tol.hor = 90/length(alpha), tol.ver = 90/length(beta), 
    cressie = FALSE, dX = numeric(0), boundaries = numeric(0), 
    cloud = FALSE, trend.beta = NULL, debug.level = 1, ...) 
{
    id1 = id2 = 0
    ret = NULL
    if (missing(cutoff)) 
        cutoff = numeric(0)
    else if (width <= 0) 
        stop("argument width should be positive")
    if (missing(width)) 
        width = numeric(0)
    if (cloud == TRUE) 
        width = 0
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
            if (!is.null(trend.beta) && length(trend.beta) > 
                0) 
                t.beta = trend.beta[[i]]
            else t.beta = numeric(0)
            .Call("gstat_new_data", as.vector(object[[i]]), 
				as.vector(locations[[i]]), as.vector(Xloc), 
				as.integer(1), as.vector(t.beta), as.integer(-1),
				as.numeric(-1), as.integer(1), numeric(0)
				, PACKAGE = "gstat"
				)
        }
    }
    else {
        nvars = 1
        if (missing(X)) 
            X = rep(1, length(object))
		if (is.null(trend.beta))
			trend.beta = numeric(0)
        .Call("gstat_new_data", as.vector(object), as.vector(locations), 
            as.vector(X), as.integer(1), as.vector(trend.beta), as.integer(-1), 
			as.numeric(-1), as.integer(1), numeric(0)
			, PACKAGE = "gstat"
			)
    }
    pos = 0
    ids = NULL
    for (id1 in nvars:1) {
        for (id2 in 1:id1) {
            if (is.null(id.names)) 
                id = ifelse(id1 == id2, paste(id1), paste(id2, 
                  ".", id1, sep = ""))
            else id = ifelse(id1 == id2, paste(id.names[id1]), 
                paste(id.names[id2], ".", id.names[id1], sep = ""))
            for (a in alpha) {
                for (b in beta) {
                  direction = as.numeric(c(a, b, tol.hor, tol.ver))
                  ret.call = .Call("gstat_variogram", 
				  		as.integer(c(id1 - 1, id2 - 1)), 
						as.numeric(cutoff), as.numeric(width), 
                    	as.numeric(direction), as.integer(cressie), 
                    	as.numeric(dX), as.numeric(boundaries)
						, PACKAGE = "gstat"
						)
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
                }
            }
        }
    }
    .Call("gstat_exit", NULL, PACKAGE = "gstat")
    ret$id = factor(ids, levels = unique(ids))
    if (cloud) 
        class(ret) = c("variogram.cloud", "data.frame")
    else class(ret) = c("variogram", "data.frame")
    ret
}
