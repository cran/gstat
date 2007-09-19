# $Id: print.variogramCloud.q,v 1.4 2007-10-18 10:13:13 edzer Exp $

"print.variogramCloud" <-
function (x, ...) 
{
	.BigInt = attr(x, ".BigInt")
	x$left = x$np %% .BigInt + 1
	x$right = x$np %/% .BigInt + 1
	x$np = NULL
    print(data.frame(x), ...)
}
