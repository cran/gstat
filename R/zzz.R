".First.lib" <-
function(lib, pkg) {
	require(lattice)
	library.dynam("gstat", pkg, lib)
	.Call("gstat_init", as.integer(1), PACKAGE = "gstat")
}
variogram <- function(object, ...) UseMethod("variogram")
