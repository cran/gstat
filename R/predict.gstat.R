"predict.gstat" <-
function (object, newdata, block = numeric(0), nsim = 0, debug.level = 1, 
    ...) 
{
    if (missing(object) || length(object$data) < 1) 
        stop("no data available")
    if (class(object) != "gstat") 
        stop("first argument should be of class gstat")
    .Call("gstat_init", as.integer(debug.level))
    nvars = length(object$data)
    pos = 1
    new.X = NULL
    names.vars = rep("xx", nvars * 2 + nvars * (nvars - 1)/2)
    for (i in 1:length(object$data)) {
        name = names(object$data)[i]
        d = object$data[[i]]
        if (d$nmax == Inf) 
            nmax = as.integer(-1)
        else nmax = as.integer(d$nmax)
        if (d$dummy) {
            tr = terms(d$locations)
            dim = length(attr(tr, "term.labels")) + attr(tr, 
                "intercept")
            .Call("gstat_new_dummy_data", as.integer(dim), as.numeric(d$beta), 
                nmax)
        }
        else {
            raw = gstat.formula(d$formula, d$locations, d$data)
            .Call("gstat_new_data", raw$y, as.vector(raw$locations), 
                as.vector(raw$X), as.numeric(d$beta), nmax)
        }
        if (!is.null(object$model[[name]])) 
            load.variogram.model(object$model[[name]], c(i - 1, i - 1))
        raw = gstat.formula.predict(d$formula, d$locations, newdata)
        if (is.null(new.X)) 
            new.X = raw$X
        else new.X = cbind(new.X, raw$X)
        names.vars[1 + (i - 1) * 2] = paste(names(object$data)[i], 
            "pred", sep = ".")
        names.vars[2 + (i - 1) * 2] = paste(names(object$data)[i], 
            "var", sep = ".")
        if (i > 1) {
            for (j in 1:(i - 1)) {
                cross = paste(names(object$data)[j], name, sep = ".")
                if (!is.null(object$model[[cross]])) 
                  load.variogram.model(object$model[[cross]], 
                    c(i - 1, j - 1))
                names.vars[nvars * 2 + pos] = paste("cov", cross, 
                  sep = ".")
                pos = pos + 1
            }
        }
    }
    if (!is.null(object$set)) 
        gstat.load.set(object$set)
    block.cols = numeric(0)
    if (!is.null(dim(block))) 
        block.cols = dim(block)[2]
    if (nsim) {
        perm = sample(1:(dim(new.X)[1]))
        ret = .Call("gstat_predict", as.integer(dim(new.X)[1]), 
            as.vector(raw$locations[perm, ]), as.vector(new.X[perm, 
                ]), as.integer(block.cols), as.vector(block), 
            as.integer(nsim))[[1]]
        ret = data.frame(cbind(raw$locations, ret[order(perm), 
            ]))
    }
    else {
        ret = .Call("gstat_predict", as.integer(dim(new.X)[1]), 
            as.vector(raw$locations), as.vector(new.X), as.integer(block.cols), 
            as.vector(block), as.integer(nsim))[[1]]
        ret = data.frame(cbind(raw$locations, ret))
    }
    .Call("gstat_exit", NULL)
    if (abs(nsim) > 0) {
        names.vars = names(object$data)
        if (length(names.vars) > 1) 
            names.vars = paste(rep(names.vars, each = abs(nsim)), 
                paste("sim", 1:abs(nsim), sep = ""), sep = ".")
        else names.vars = paste("sim", 1:abs(nsim), sep = "")
        if (abs(nsim) == 1) 
            ret = ret[, 1:(dim(ret)[2] - 1)]
    }
    names(ret) = c(dimnames(raw$locations)[[2]], names.vars)
    return(ret)
}
