"get.direction.unitv" <- function(alpha, beta) {
	cb = cos(beta)
	c(cb * sin(alpha), cb * cos(alpha), sin(beta))
}

"vgm.panel.xyplot" <-
function (x, y, type = "p", pch = plot.symbol$pch, col, 
	col.line = plot.line$col, col.symbol = plot.symbol$col, lty = 
	plot.line$lty, cex = plot.symbol$cex, lwd = plot.line$lwd, 
	model = model, direction = direction, labels, shift = shift, ...) 
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (length(x) > 0) {
        if (!missing(col)) {
            if (missing(col.line)) 
                col.line <- col
            if (missing(col.symbol)) 
                col.symbol <- col
        }
        plot.symbol <- trellis.par.get("plot.symbol")
        plot.line <- trellis.par.get("plot.line")
        lpoints(x = x, y = y, cex = cex, col = col.symbol, pch = pch)
        if (!missing(model) && !is.null(model)) {
            ang.hor <- pi * (direction[1]/180)
			ang.ver <- pi * (direction[2]/180)
            dir <- get.direction.unitv(ang.hor, ang.ver)
            ret <- variogram.line(model, max(x), dir = dir)
            llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, 
                lwd = lwd)
        }
		if (!is.null(labels))
			ltext(x = x + shift * max(x), y = y, labels = labels)
    }
}

"xvgm.panel.xyplot" <-
function (x, y, subscripts, type = "p", pch = plot.symbol$pch, 
    col, col.line = plot.line$col, col.symbol = plot.symbol$col, 
    lty = plot.line$lty, cex = plot.symbol$cex, ids, lwd = plot.line$lwd, 
    model = model, direction = direction, labels, shift = shift, ...) 
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    id <- as.character(ids[subscripts][1])
    if (length(x) > 0) {
        if (!missing(col)) {
            if (missing(col.line)) 
                col.line <- col
            if (missing(col.symbol)) 
                col.symbol <- col
        }
        plot.symbol <- trellis.par.get("plot.symbol")
        plot.line <- trellis.par.get("plot.line")
        lpoints(x = x, y = y, cex = cex, col = col.symbol, pch = pch)
        if (!missing(model) && !is.null(model)) {
			if (inherits(model, "gstat"))
				m = model$model
			else
				m = model
			if (!is.list(m))
				stop("model argument not of class gstat or list")
			if (is.list(m) && !is.null(m[[id]])) {
                ang.hor <- pi * (direction[1]/180)
				ang.ver <- pi * (direction[2]/180)
                dir <- get.direction.unitv(ang.hor, ang.ver)
				ret <- variogram.line(m[[id]], max(x), dir = dir)
				llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, 
					lwd = lwd)
			}
        }
        if (!is.null(labels)) 
            ltext(x = x + shift * max(x), y = y, labels = labels[subscripts])
    }
}

"vgm.dir.panel.xyplot" <-
function (x, y, subscripts, type = "p", pch = plot.symbol$pch, 
    col, col.line = plot.line$col, col.symbol = plot.symbol$col, 
    lty = plot.line$lty, cex = plot.symbol$cex, lwd = plot.line$lwd, 
    model = model, dir.hor = dir.hor, labels, shift = shift, ...) 
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (length(x) > 0) {
        if (!missing(col)) {
            if (missing(col.line)) 
                col.line <- col
            if (missing(col.symbol)) 
                col.symbol <- col
        }
        plot.symbol <- trellis.par.get("plot.symbol")
        plot.line <- trellis.par.get("plot.line")
        lpoints(x = x, y = y, cex = cex, col = col.symbol, pch = pch)
        if (!missing(model) && !is.null(model)) {
            dir <- c(1, 0, 0)
            if (!missing(dir.hor)) {
                ang.hor <- pi * (dir.hor[subscripts][1]/180.0)
                dir <- get.direction.unitv(ang.hor, 0)
            }
            ret <- variogram.line(model, max = max(x), dir = dir)
            llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, 
                lwd = lwd)
        }
		if(!is.null(labels))
			ltext(x + shift * max(x), y, labels = labels[subscripts])
    }
}
