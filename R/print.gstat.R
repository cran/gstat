"print.gstat" <-
function (x, ...) 
{
    if (missing(x) || !inherits(x, "gstat"))
        stop("wrong call")
    data.names <- names(x$data)
    if (length(data.names)) 
        cat("data:\n")
    for (n in data.names) {
        fstr = paste(x$data[[n]]$formula[c(2, 1, 3)], collapse = "")
        lstr = paste(x$data[[n]]$locations[c(1, 2)], collapse = "")
        cat(n, ": formula =", fstr, "; locations =", lstr, ";")
        if (!is.null(x$data[[n]]$data)) {
            data.dim = dim(x$data[[n]]$data)
            cat(" data dim =", data.dim[1], "x", data.dim[2])
        }
        else {
            if (x$data[[n]]$dummy) 
                cat(" dummy data")
            else cat(" NULL data")
        }
        cat("\n")
    }
    xx.names = xx = NULL
    for (n in data.names) {
        m = x$model[[n]]
        if (!is.null(m)) {
            xx = rbind(xx, m)
            if (nrow(m) == 1) 
                xx.names = c(xx.names, n)
            else xx.names = c(xx.names, paste(n, "[", 1:nrow(m), 
                "]", sep = ""))
        }
    }
    if (length(data.names) > 1) {
        for (j in 2:length(data.names)) {
            for (i in 1:(j - 1)) {
                n = paste(data.names[i], data.names[j], sep = ".")
                m = x$model[[n]]
                if (!is.null(m)) {
                  xx = rbind(xx, m)
                  if (nrow(m) == 1) 
                    xx.names = c(xx.names, n)
                  else xx.names = c(xx.names, paste(n, "[", 1:nrow(m), 
                    "]", sep = ""))
                }
            }
        }
    }
    if (!is.null(xx)) {
        cat("variograms:\n")
        row.names(xx) = xx.names
        print(xx, ...)
    }
    if (!is.null(x$set)) {
        s = gstat.set(x$set)
        for (i in 1:length(s)) cat(s[i], "\n")
    }
    invisible(x)
}
