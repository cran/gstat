.onLoad <-
function(lib, pkg) {
	library.dynam("gstat", pkg, lib)
	.Call("gstat_init", as.integer(1), PACKAGE = "gstat")
	require(lattice)
}


variogram <- function(object, ...) UseMethod("variogram")
