"gstat.formula.predict" <-
function (formula, locations, data) 
{
    delete.resp.f <- function(x) {
        f <- formula(x)
        if (length(f) == 3) 
            f[[2]] <- NULL
        f
    }
    call = match.call()
    X = model.matrix(delete.resp.f(formula), data = data)
    if (nrow(X) == 1 && nrow(data) > 1) 
        X = as.matrix(rep(X, nrow(data)))
    m = match.call(expand = FALSE)
    m$method = m$model = m$x = m$y = m$... = NULL
    m$formula = m$locations
    m$locations = NULL
    m[[1]] = as.name("model.frame")
    m = eval(m, sys.frame(sys.parent()))
    Terms = attr(m, "terms")
    attr(Terms, "intercept") = 0
    # locations = model.matrix(locations, data = data)
    locations = model.matrix(Terms, m)
    list(locations = locations, X = X, call = call)
}
