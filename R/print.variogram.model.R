"print.variogram.model" <-
function (x, ...) 
{
    df <- data.frame(x)
	if (!any(df[, "model"] == "Mat"))
		df$kappa <- NULL
    if (!any(df[, "anis2"] != 1))  {
		df$anis2 <- NULL
		df$ang2 <- NULL
		df$ang3 <- NULL
		if (!any(df[, "anis1"] != 1))  {
			df$anis1 <- NULL
			df$ang1 <- NULL
		}
	} 
    print(df, ...)
	invisible(x)
}
