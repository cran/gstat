"xyz2img" <-
function (xyz, zcol = 3, xcol = 1, ycol = 2) 
{
    if (ncol(xyz) < 3) 
        stop("xyz object should have at least three columns")
    z = xyz[, zcol]
    x = xyz[, xcol]
    y = xyz[, ycol]
    xx = sort(unique(x))
    yy = sort(unique(y))
    nx = length(xx)
    ny = length(yy)
    nmax = max(nx, ny)
    difx = diff(xx)
    if (diff(range(unique(difx))) > 1e-15) 
        stop("x intervals are not constant")
    dify = diff(yy)
    if (diff(range(unique(dify))) > 1e-15) 
        stop("y intervals are not constant")
    dx = difx[1]
    dy = dify[1]
    ratio = (nx * dx)/(ny * dy)
    xmin = min(xx)
    xmax = max(xx)
    xrange = xmax - xmin
    ymin = min(yy)
    ymax = max(yy)
    yrange = ymax - ymin
    zz = matrix(NA, nrow = nx, ncol = ny)
    xx = seq(xmin, xmax, dx)
    yy = seq(ymin, ymax, dy)
    row = round((x - xmin)/dx) + 1
    col = round((y - ymin)/dy) + 1
    for (i in 1:length(x)) zz[row[i], col[i]] = z[i]
    list(x = xx, y = yy, z = zz)
}
