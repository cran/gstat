gstat.set <- function(set) {
	if(!is.list(set))
		stop("set should be a list")
	if (length(set) == 0)
		return(NULL)
	ret = NULL
	n = names(set)
	for (i in (1:length(set))) {
		val = set[[i]]
		if (is.character(val))
			val = paste("'", val, "'", sep = "")
		str = paste("set ", n[i], " = ", val, ";", sep="")
		ret = c(ret, str)
	}
	ret
}

gstat.load.set <- function(set) {
	str = gstat.set(set)
	if (!is.null(str)) {
		ret = .C("Cload_gstat_command", str, as.integer(length(str)), 
			as.integer(0))[[3]]
		if (ret != 0)
			stop(paste("error occured when parsing command:", str[ret]))
	}
	invisible()
}
