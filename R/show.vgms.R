"show.vgms" <-
function(min = 1e-12 * max, max = 3, n = 50, sill = 1, range = 1,
	models = as.character(vgm()[c(1:16)]), nugget = 0, kappa.range = 0.5,
	plot = TRUE) 
{

	zero.range.models = c("Nug", "Int", "Lin", "Err")
	# print(models)
	i = 0
	if (length(kappa.range) > 1) { # loop over kappa values for Matern model:
		data = matrix(NA, n * length(kappa.range), 2)
		v.level = rep("", n * length(kappa.range))
		for (kappa in kappa.range) {
			v = vgm(sill, "Mat", range, nugget = nugget, kappa = kappa)
			x = variogram.line(v, 0, 1, 0)
			data[(i*n+1), ] = as.matrix(x)
			x = variogram.line(v, max, n - 1, min)
			data[(i*n+2):((i+1)*n), ] = as.matrix(x)
			m.name = paste("vgm(", sill, ",\"Mat\",", range, sep = "")
			if (nugget > 0)
				m.name = paste(m.name, ",nugget=", nugget, sep = "")
			m.name = paste(m.name, ",kappa=", kappa, ")", sep = "")
			v.level[(i*n+1):((i+1)*n)] = rep(m.name, n)
			i =  i + 1
		}
	} else {
		data = matrix(NA, n * length(models), 2)
		v.level = rep("", n * length(models))
		for (m in models) {
			this.range = ifelse(!is.na(pmatch(m,zero.range.models)), 0, range)
			v = vgm(sill, m, this.range, nugget = nugget, kappa = kappa.range)
			x = variogram.line(v, 0, 1, 0)
			data[(i*n+1), ] = as.matrix(x)
			x = variogram.line(v, max, n - 1, min)
			data[(i*n+2):((i+1)*n), ] = as.matrix(x)
			m.name = paste("vgm(", sill, ",\"", m, "\",", this.range, sep = "")
			if (nugget > 0)
				m.name = paste(m.name, ",nugget=", nugget, sep = "")
			m.name = paste(m.name, ")", sep = "")
			v.level[(i*n+1):((i+1)*n)] = rep(m.name, n)
			i =  i + 1
		}
	}
	dframe = data.frame(semivariance = data[,2], distance = data[,1], 
		model = factor(v.level, levels = unique(v.level)))
	vgm.panel = function(x,y) {
		n = length(x)
		lpoints(x[1],y[1])
		llines(x[2:n],y[2:n])
	}
	if (!plot)
		dframe
	else
		xyplot(semivariance ~ distance | model, dframe, 
			panel = vgm.panel, as.table = TRUE)
}
