"image.data.frame" <-
function (x, zcol = 3, xcol = 1, ycol = 2, ...)
{
    image.default(xyz2img(xyz = x, zcol = zcol, xcol = xcol, ycol = ycol),...)
}
