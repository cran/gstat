"map.to.lev" <-
function (data, xcol = 1, ycol = 2, zcol = c(3, 4), ns = names(data)[zcol]) 
{
	# print(ns)
	len = dim(data)[1]
    d = matrix(nrow = len * length(zcol), ncol = 3)
	xnames = NULL
    if (length(ns) > 1 && length(ns) != length(zcol)) 
        stop("names should have length 1 or equal to length of zcol")
	nr = 1
    for (i in zcol) {
        if (length(ns) == 1) 
            nm = rep(paste(ns, nr), len)
        else nm = rep(ns[nr], len)
        range = (1 + (nr - 1) * len):(nr * len)
        d[range,] = cbind(data[, xcol], data[,ycol], data[, i]) 
		xnames = c(xnames, nm)
        nr = nr + 1
    }
	d = data.frame(d, xnames)
    # d$name = as.factor(d$name)
	names(d) = c("x", "y", "z", "name")
    d
}
