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
    pin = par("pin")
    pinratio = pin[1]/pin[2]
    if (ratio > pinratio) {
        yrange = yrange * ratio/pinratio
        ny.new = max(round(yrange/dy), ny)
        ymin = min(yy) - dy * round((ny.new - ny)/2)
        ymax = ymin + dy * (ny.new - 1)
        ny = ny.new
    }
    else {
        xrange = xrange * pinratio/ratio
        nx.new = max(round(xrange/dx), nx)
        xmin = min(xx) - dx * round((nx.new - nx)/2)
        xmax = xmin + dx * (nx.new - 1)
        nx = nx.new
    }
    zz = matrix(NA, nrow = nx, ncol = ny)
    xx = seq(xmin, xmax, dx)
    yy = seq(ymin, ymax, dy)
    row = round((x - xmin)/dx) + 1
    col = round((y - ymin)/dy) + 1
    for (i in 1:length(x)) zz[row[i], col[i]] = z[i]
    list(x = xx, y = yy, z = zz)
}
