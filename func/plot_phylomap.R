### FUNCTION
# plot_phylomap
### ARGUMENTS
# coords: data frame with individual identifier and geographical coordinates (whose names are partial matches of 'lat', 'lon' or 'x', 'y')
# class: data frame with individual identifier and classification into groups (species, populations, ...)
# ind: name of individual identifier (e.g. 'ID')
# group: name of group identifier (e.g. 'Species')
# color: named vector of colors, names must correspond to group labels; if missing a default palette is used
# border: color of point borders; if 'auto' it is identical to point fillings
# size: size of points specified as a percentage of the lesser of latitude / longitude ranges
# map: map layer to serve as background instead of country broders from 'maps' package
# device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
# file: name of pdf file if device == "pdf"
# width: width of graphical device
# height: height of graphical device
# mai: size of outer margins in inches, recycled if necessary
# xlim, ylim: limits of x and y axes
# axes: whether to display axes
# cex.axis: magnification to be used for axis annotation
# expr: an expression (or a list of them) allowing to add other elements (e.g. legend)

plot_phylomap <- function(coords, class, ind, group, color, border="black", size=4, lwd=1, map=NULL, device, width=7, height=7, file=NULL, mai=0.02, xlim, ylim, axes=TRUE, cex.axis=1.25, expr) {
	
	partCircle <- function(x, y, prop, size, res, col, border, lwd) {
		rs <- seq(0, 2 * pi, len=res)
		pts <- cbind(0.5 * cos(rs), 0.5 * sin(rs))
		pts <- size * pts + rep(1, res) %*% t(c(x,y))
		prop <- prop / sum(prop)
		parts <- round(quantile(seq(nrow(pts)), probs=c(0,cumsum(prop))))
		if (is.null(names(border)) | !all(names(border) %in% names(col))) {
			border <- setNames(rep_len(border, length(col)), names(col))
		}
		if (all(names(prop) %in% names(col)) & !is.null(names(prop))) {
			col <- col[names(prop)]
			border <- border[names(prop)]
		} else {
			col <- rep_len(col, length(prop))
			border <- rep_len(border, length(prop))
		}
		for (i in seq_along(prop)) {
			ipts <- rbind(c(x,y), pts[parts[i]:parts[i+1],], c(x,y))
			polygon(ipts, border="transparent", col=col[i], lwd=1)
			polygon(ipts[-c(1,nrow(ipts)),], border=border[i], col="transparent", lwd=lwd)
		}
	}
	
	if (missing(ind)) {
		ind <- ifelse(names(class)[1] %in% names(coords), names(class)[1], names(coords)[1])
	}
	if (missing(group)) {
		group <- setdiff(names(class),ind)[1]
	}
	xy <- pmatch(c("lon","lat"), tolower(names(info)))
	if (any(is.na(xy))) {
		xy <- pmatch(c("x","y"), tolower(names(info)))
	}
	
	coords <- coords[coords[,ind] %in% class[,ind],]
	m <- match(coords[,ind], class[,ind])
	coords$OTU <- class[m, group]	
	xy <- coords[,xy]
	sp <- coords[,"OTU"]
	L <- sort(unique(sp))
	K <- length(L)
	if (missing(color)) {
		palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")
		color <- setNames(palette[seq(K)], L)
	}

	if (missing(xlim)) {
		xlim <- range(xy[,1])
	}
	if (missing(ylim)) {
		ylim <- range(xy[,2])
	}
	circ <-  0.04 * min(diff(xlim), diff(ylim)) * 0.5 * cos(0)
	xlim <- xlim + c(-1, 1) * circ
	ylim <- ylim + c(-1, 1) * circ
	loc <- unique(xy)
	sites <- vector("list", nrow(loc))
	for (i in seq_along(sites)) {
		sites[[i]] <- which(xy[,1] == loc[i,1] & xy[,2] == loc[i,2])
	}
	prop <- lapply(sites, function(x, sp) table(sp[x]), sp=sp)
	ord <- order(sapply(prop, length), sapply(prop, sum))
	loc <- loc[ord,]
	prop <- prop[ord]
	
	# plotting
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "phylo_map.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}
	mai <- rep_len(mai, 4)	
	if (isTRUE(axes)) {
		onechar <- strheight(s="L", units="inches", cex=1) * cex.axis
		mai[1] <- 3 * onechar + 0.5 * onechar / 0.2
		mai[2] <- 3 * onechar + 0.5 * onechar / 0.2
	}
	
	par(mai=mai)
	plot(xlim, ylim, asp=1, type="n", bty="n", axes=FALSE, ann=FALSE)
	if (!is.null(map)) {
		plot(map, asp=1, xlim=xlim, ylim=ylim, add=TRUE)
	} else {
		maps::map(xlim=xlim, ylim=ylim, add=TRUE)
	}

	size <- 0.01 * size * min(diff(xlim), diff(ylim))
	if (border == "auto") {
		border <- lapply(prop, function(x) color[names(x)])
	} else {
		border <- rep(border, length(prop))
	}
	for (i in seq_along(sites)) {
		partCircle(x=loc[i,1], y=loc[i,2], prop=prop[[i]], size=size, res=200, col=color[names(prop[[i]])], border=border[[i]], lwd=lwd)
	}

	if (isTRUE(axes)) {
		line <- 0.5 * onechar / 0.2
		axis(side=1, cex.axis=cex.axis, line=line)
		axis(side=2, cex.axis=cex.axis, line=line)
		title(xlab="Longitude", ylab="Latitude", cex.lab=cex.axis)
	}

	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		} else if (all(sapply(expr, is.expression))) {
			for (i in seq_along(expr)) {
				eval(expr[[i]])
			}
		}		
	}

# closing pdf device	
	if (isTRUE(device == "pdf")) {
		dev.off()
	}

}
