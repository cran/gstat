"bubble" <-
function(data, xcol = 1, ycol = 2, zcol = 3, fill = TRUE, maxsize = 3,
		do.sqrt = TRUE, pch, col = c(2,3), ...) {
	x = data[,xcol]
	y = data[,ycol]
	z = data[,zcol]
	d = data.frame(x=x, y=y)
	if (missing(pch))
		pch = ifelse(fill, 16, 1)
	col = ifelse(z < 0, col[1], col[2])
	az = abs(z)
	if (do.sqrt)
		az = sqrt(az)
	cex = maxsize * az/max(az)
	xyplot(y~x, d, col = col, cex = cex, pch = pch, asp = mapasp(d), ...)
}
