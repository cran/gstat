"vgm.dir.panel.xyplot" <-
function (x, y, subscripts, type = "p", pch = plot.symbol$pch, 
    col, col.line = plot.line$col, col.symbol = plot.symbol$col, 
    lty = plot.line$lty, cex = plot.symbol$cex, lwd = plot.line$lwd, 
    model = model, dir.hor = dir.hor, labels, ...) 
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
        if (!missing(model)) {
            if (!missing(dir.hor)) {
                ang <- dir.hor[subscripts]
                ang.pi <- pi * (ang[1]/180)
                dir <- c(sin(ang.pi), cos(ang.pi), 0)
            }
            else dir <- c(1, 0, 0)
            ret <- variogram.line(model, max = max(x), dir = dir)
            llines(x = ret$dist, y = ret$gamma, lty = lty, col = col.line, 
                lwd = lwd)
        }
		if(!is.null(labels))
			ltext(x + 0.05 * max(x), y, labels = labels)
    }
}
