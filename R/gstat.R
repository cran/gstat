"gstat" <-
function (g, id, formula, locations, data = NULL, model = NULL, 
    beta, nmax = Inf, dummy = FALSE, set, fill.all = FALSE) 
{
	if (fill.all) {
		if (missing(g) || is.null(model))
			stop("fill.all assumes object g and model are supplied")
        g.names = names(g$data)
		for (i in 1:length(g.names)) {
           	g$model[[paste(g.names[i])]] = model
			for (j in (i+1):length(g.names))
           		g$model[[paste(g.names[i], g.names[j], sep = ".")]] = model
		}
        return(g)
	} 
    if (!missing(g) && class(g) == "gstat" && !missing(id) && 
        !missing(model) && missing(formula) && missing(locations)) {
        g.names = names(g$data)
		if (length(id) == 2) {
           	m1 = match(id[1], g.names)
           	m2 = match(id[2], g.names)
           	if (is.na(m1)) 
               	stop("first id does not match available data")
           	if (is.na(m1)) 
               	stop("second id does not match available data")
           	nm = paste(g.names[min(m1, m2)], g.names[max(m1, 
               	m2)], sep = ".")
        }
        else if (is.na(match(id, g.names))) 
           	stop("id does not match available data")
        g$model[[nm]] = model
        return(g)
    }
    if (class(formula) != "formula") 
        stop("argument formula should be of class formula")
    if (class(locations) != "formula") 
        stop("argument locations should be of class formula")
    if (missing(beta) || is.null(beta)) 
        beta = numeric(0)
    if (missing(g)) {
        g = list()
        g[["data"]] = list()
        g[["model"]] = list()
    }
    if (missing(id)) 
        id = paste("var", length(g$data) + 1, sep = "")
    g$data[[id]] = list(formula = formula, locations = locations, 
        data = data, beta = beta, nmax = nmax, dummy = dummy)
    g$model[[id]] = model
    if (!missing(set)) {
        if (!is.list(set)) 
            stop("argument set should be a list")
        g$set = set
    }
    class(g) = c("gstat", "list")
    g
}
