"print.variogram.model" <-
function (x, ...) 
{
    df <- data.frame(x)
    if (any(df[, "anis2"] != 1)) 
        print(df, ...)
    else if (any(df[, "anis1"] != 1)) 
		print(df[, c(1:3,4,7)], ...)
	else print(df[, 1:3], ...)
}
