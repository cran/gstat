"plot.variogram" <-
function (x, ylim, xlim, xlab = "distance", ylab = "semivariance", 
    multipanel = TRUE, plot.numbers = FALSE, scales, ids = x$id, 
    ...) 
{
    if (missing(ylim) && min(x$gamma) >= 0) 
        ylim = c(0, 1.04 * max(x$gamma))
    if (missing(xlim)) 
        xlim = c(0, 1.04 * max(x$dist))
    labels = NULL
    if (plot.numbers == TRUE) 
        labels = as.character(x$np)
    if (length(unique(x$dir.hor)) > 1) {
        if (multipanel) 
            xyplot(gamma ~ dist | as.factor(dir.hor), subscripts = TRUE, 
                panel = vgm.dir.panel.xyplot, data = x, xlim = xlim, 
                ylim = ylim, xlab = xlab, ylab = ylab, dir.hor = x$dir.hor, 
                labels = labels, ...)
        else {
            pch = as.integer(as.factor(x$dir.hor))
            xyplot(gamma ~ dist, data = x, type = c("p", "l"), 
                groups = pch, xlim = xlim, ylim = ylim, xlab = xlab, 
                ylab = ylab, pch = pch, ...)
        }
    }
    else if (length(unique(ids)) > 1) {
        n = floor(sqrt(2 * length(unique(x$id))))
        skip = NULL
        for (row in n:1) for (col in 1:n) skip = c(skip, row < 
            col)
        if (missing(scales)) 
            scales = list(y = list(relation = "free"))
        xyplot(gamma ~ dist | as.factor(id), data = x, xlim = xlim, 
            ylim = ylim, xlab = xlab, ylab = ylab, ids = ids, 
            panel = xvgm.panel.xyplot, labels = labels, scales = scales, 
            layout = c(n, n), skip = skip, prepanel = function(x, 
                y) list(ylim = c(min(0, y), max(0, y))), ...)
    }
    else xyplot(gamma ~ dist, data = x, panel = vgm.panel.xyplot, 
        xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, labels = labels, 
        ...)
}
