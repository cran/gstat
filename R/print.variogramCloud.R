# $Id: print.variogramCloud.q,v 1.2 2006-02-10 19:01:07 edzer Exp $

"print.variogramCloud" <-
function (x, ...) 
{
	x$left = x$np %% 2^16 + 1
	x$right = x$np %/% 2^16 + 1
	x$np = NULL
    print(data.frame(x), ...)
}
