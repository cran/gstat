"bubble" <-
function (data, xcol = 1, ycol = 2, zcol = 3, fill = TRUE, maxsize = 3, 
    do.sqrt = TRUE, pch, col = c(2, 3), key.entries = quantile(data[,zcol]),
    ...) 
{
    x = data[, xcol]
    y = data[, ycol]
    z = data[, zcol]
    d = data.frame(x = x, y = y)
    if (missing(pch)) 
        pch = ifelse(fill, 16, 1)
    z.col = ifelse(z < 0, col[1], col[2])
    q = key.entries
    q.pch = rep(pch, length(q))
    q.text = as.character(round(q, 3))
    q.col = ifelse(q < 0, col[1], col[2])
    az = abs(z)
    q = abs(q)
    if (do.sqrt) {
        az = sqrt(az)
	q = sqrt(q)
    }
    cex = maxsize * az/max(az)
    q.cex = maxsize * q/max(az)

    key = list(space = "right", points = list(pch = q.pch, col = q.col, 
    	cex = q.cex), text = list(q.text))
    xyplot(y ~ x, d, col = z.col, cex = cex, pch = pch, asp = mapasp(d), 
        key = key, ...)
}
