"bpy.colors" <-
function (n = 100) 
{
    if ((n <- as.integer(n[1])) > 0) {
        i = (0:(n - 1))/(n - 1)
		r = ifelse(i < .25, 0, ifelse(i < .57, i/.32 - .78125, 1))
		g = ifelse(i < .42, 0, ifelse(i < .92, 2*i - .84, 1))
		b = ifelse(i < .25, 4*i, ifelse(i < .42, 1, 
			ifelse(i < .92, -2 * i + 1.84, i/.08 - 11.5)))
        rgb(r, g, b)
    }
    else character(0)
}
