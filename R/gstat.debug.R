"gstat.debug" <-
function(level = 0)
{
	invisible(.Call("gstat_debug_level", as.integer(level)))
}
