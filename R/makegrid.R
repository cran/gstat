makegrid = function(x, y, n = 10000, nsig = 2, margin = 1.05, cell.size) {
        dx = range(x)
        dy = range(y)
		if (missing(cell.size)) {
        	cell.area = margin * diff(dx) * margin * diff(dy) / n
        	cell.size = signif(sqrt(cell.area), nsig)
		} 
        midx = signif(mean(dx), nsig)
        midy = signif(mean(dy), nsig)
        nx = ceiling((midx - min(x))/cell.size) * 2 + 1
        ny = ceiling((midy - min(y))/cell.size) * 2 + 1
        minx = midx - ceiling((midx - min(x))/cell.size) * cell.size
        maxx = midx + ceiling((midx - min(x))/cell.size) * cell.size
        miny = midy - ceiling((midy - min(y))/cell.size) * cell.size
        maxy = midy + ceiling((midy - min(y))/cell.size) * cell.size
        expand.grid(x = seq(minx, maxx, length = nx),
			y = seq(miny, maxy, length = ny))
}
