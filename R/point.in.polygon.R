"point.in.polygon" <-
function(point.x, point.y, pol.x, pol.y) {
	as.logical(.Call("gstat_pip", 
		as.numeric(point.x),
		as.numeric(point.y),
		as.numeric(pol.x),
		as.numeric(pol.y)
		, PACKAGE = "gstat")
		)
}
