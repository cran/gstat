"gstat.formula" <-
function (formula, locations, data) 
{
    call = match.call()
    m = match.call(expand = FALSE)
    Y = NULL
    m$method = m$model = m$x = m$y = m$... = NULL
    m$drop.unused.levels = TRUE
    m[[1]] = as.name("model.frame")
    m$locations = NULL
    m = eval(m, sys.frame(sys.parent()))
    Terms = attr(m, "terms")
    Y = model.extract(m, response)
    if (length(Y) == 0) 
        stop("no response variable present in formula")
    X = model.matrix(Terms, m)
    if (class(locations) != "formula") 
        stop("locations argument should be a formula, such as ~x+y")
    m = match.call(expand = FALSE)
    m$method = m$model = m$x = m$y = m$... = NULL
    m$formula = locations
    m$locations = NULL
    m[[1]] = as.name("model.frame")
    m = eval(m, sys.frame(sys.parent()))
    Terms = attr(m, "terms")
    attr(Terms, "intercept") = 0
    if ((yvar = attr(Terms, "response")) > 0) 
        stop("no response allowed in locations formula")
    locations = model.matrix(Terms, m)
    list(y = Y, locations = locations, X = X, call = call)
}
