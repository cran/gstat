".First.lib" <-
function(lib, pkg)
{
	require(lattice)
	library.dynam("gstat", pkg, lib)
	.Call("gstat_init", as.integer(1))
}

variogram <- function(object, ...) UseMethod("variogram")

# ".on.attach" <-
# function() {
# 	.Call("gstat_init")
# 	print("gstat initialized")
# }
