"print.variogram.cloud" <-
function (x, ...) 
{
	x$left = x$np %% 2^16 + 1
	x$right = x$np %/% 2^16 + 1
	x$np = NULL
    print(data.frame(x), ...)
}
