"fit.lmc" <-
function (v, g, model, fit.ranges = FALSE, fit.lmc = !fit.ranges, 
    ...) 
{
    posdef = function(X) {
        q = eigen(X)
        d = q$values
        d[d < 0] = 0
        q$vectors %*% diag(d) %*% t(q$vectors)
    }
    if (!inherits(v, "variogram"))
        stop("v should be of class variogram")
    if (!inherits(g, "gstat"))
        stop("g should be of class gstat")
    if (!missing(model)) {
        if (!inherits(model, "variogram.model"))
            stop("model should be of class variogram.model")
    }
    n = names(g$data)
    for (i in 1:length(n)) {
        for (j in i:length(n)) {
            name = ifelse(i == j, n[i], paste(n[i], n[j], sep = "."))
            x = v[v$id == name, ]
            if (nrow(x) == 0) 
                stop(paste("variogram", name, "not present"))
            m = g$model[[name]]
            if (!missing(model)) 
                m = model
            g$model[[name]] = fit.variogram(x, m, fit.ranges = fit.ranges, 
                ...)
        }
    }
    if (fit.lmc) {
        m = g$model[[n[1]]]
        for (k in 1:nrow(m)) {
            psill = matrix(NA, nrow = length(n), ncol = length(n))
            for (i in 1:length(n)) {
                for (j in i:length(n)) {
                  name = ifelse(i == j, n[i], paste(n[i], n[j], 
                    sep = "."))
                  psill[i, j] = psill[j, i] = g$model[[name]][k, 
                    "psill"]
                }
            }
            psill = posdef(psill)
            for (i in 1:length(n)) {
                for (j in i:length(n)) {
                  name = ifelse(i == j, n[i], paste(n[i], n[j], 
                    sep = "."))
                  g$model[[name]][k, "psill"] = psill[i, j]
                }
            }
        }
    }
    g
}
