### FUNCTION
# plot_clades
### ARGUMENTS
# phy: tree object of class 'phylo'
# clades: Vector of clade labels (order must correspond to 'phy$tip.label') or data frame (or matrix) with tip labels in the first column and clade assignment in the second. The classification need not to be comprehensive.
# color: named vector of colors. Names correspond to either clade labels or tip labels (if coloration is unlinked from clades)
# stem: color of stem (unclassified) branches
# type: tree type, either phylogram (default) or cladogram
# lines: logical, whether to display vertical lines
# order: either logical (whether to alphabetically order clades) or character (clade labels in the required order)
# edge.color: branch colors, vector of the same length as nrow(phy$edge), if specifies it supersedes color in painting of branches, while color of clade lines and names is not affected
# edge.width: branch width
# scale: whether to display scale bar as x-axis
# tip_labels: named vector of tip labels to be displayed (names indicate tip numbers)
# node_labels: named vector of node labels to be displayed (names indicate node numbers, i.e. higher than Ntip(phy))
# edge_labels: named vector of edge labels to be displayed (names indicate edge order in phy$edge)
# root_label: root label to be added to node_labels (e.g., "100" when indicating bootstrap support)
# anc_labels: logical, whether to show node_labels only for clade ancestors
# device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
# file: name of pdf file if device == "pdf"
# width: width of graphical device
# height: height of graphical device
# mai: size of outer margins in inches, recycled if necessary
# extension: extension of axes expressed as percents of the data range
# line.pos: position of clade lines (negative values ~ shift to the left)
# line.over: overhang of clade lines (as a proportion of single character height)
# name.col: color of clade names
# name.cex: size of clade names
# name.font: font of clade names
# name.las: orientation of clades names relative to the axis (default: perpendicular)
# tip.lab.pos: gap between tree tips and tip labels as a proportion of the maximum root-tip distance
# tip.lab.cex: size of tip labels
# tip.lab.font: font of tip labels
# node.lab.pos: node label position as a proportion of the supporting branches
# node.lab.cex: size of node label point
# node.lab.lwd: width of node label point border
# node.lab.asp: shape of node label, aspect ratio of an ellipse, default=1 means cricle
# node.lab.col: color of node label point border
# node.lab.bg: color of node label point background
# edge.lab.pos: edge label position as a proportion of the branch from its start (= root side)
# edge.lab.cex: size of edge label point
# edge.lab.lwd: width of edge label point border
# edge.lab.asp: shape of edge label, aspect ratio of an ellipse, default=1 means cricle
# edge.lab.col: color of edge label point border
# edge.lab.bg: color of edge label point background
# text.lab.pos: adjustment of text position in node / edge label (adj of 'text')
# text.lab.cex: size of text in node / edge label
# text.lab.font: font of text in node / edge label
# text.lab.col: color of text in node / edge label

plot_clades <- function(phy, clades, color, stem="black", type="phylogram", lines=TRUE, names=TRUE, order=FALSE, edge.color, edge.width=3, scale=TRUE, tip_labels=TRUE, node_labels=TRUE, edge_labels=NULL, root_label=NULL, anc_labels=TRUE, device, width=7, height=7, file=NULL, mai, extension=8, line.pos=0, line.over=0.5, name.col=NULL, name.cex=1.25, name.font=1, name.las=2, tip.lab.pos=0.005, tip.lab.cex=1, tip.lab.font=3, node.lab.pos=1, node.lab.cex=1, node.lab.lwd=1, node.lab.asp=1, node.lab.col="black", node.lab.bg="white", edge.lab.pos=0.5, edge.lab.cex=1, edge.lab.lwd=1, edge.lab.asp=1, edge.lab.col="black", edge.lab.bg="white", text.lab.pos=c(0.5,0.5), text.lab.cex=1, text.lab.font=1, text.lab.col="black") {

	# functions
	getPaths <- function(from, to, phy) lapply(to, function(x) ape::nodepath(phy, from=from, to=x))
	getOffspring <- function(tip, phy) match(tip, phy$tip.label)
	getAncestor <- function(tip, phy) {
		ifelse(length(tip) == 1, phy$edge[phy$edge[,2] == match(tip, phy$tip.label), 1], ape::getMRCA(phy, tip))
	}
	makeEdges <- function(paths) {
		unique(do.call(rbind, lapply(paths, function(x) cbind(x[-length(x)], x[-1]))))
	}
	drawPhylogram <- function(phy, xy, lwd=1, col=1, ...) {
		lwd <- rep_len(lwd, nrow(phy$edge))
		col <- rep_len(col, nrow(phy$edge))
		for (i in seq(nrow(phy$edge))) {
			edg <- xy[phy$edge[i,],]
			lines(matrix(edg[c(1,2,4,4)],2,2), lwd=lwd[i], col=col[i], ...)
			lines(matrix(edg[c(1,1,3,4)],2,2), lwd=lwd[i], col=col[i], ...)
		}
	}
	drawCladogram <- function(phy, xy, lwd=1, col=1, ...) {
		lwd <- rep_len(lwd, nrow(phy$edge))
		col <- rep_len(col, nrow(phy$edge))
		for (i in seq(nrow(phy$edge))) {
			lines(xy[phy$edge[i,],], lwd=lwd[i], col=col[i], ...)
		}
	}
	drawTree <- list(phylogram=drawPhylogram, cladogram=drawCladogram)[[type]]
	drawEllipse <- function(x, y, size, w, h, r=0, res=200, col, border, lwd) {
		rs <- seq(0, 2 * pi, len=res)
		pts <- cbind(0.5 * w * cos(rs), 0.5 * h * sin(rs))
		rot <- matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
		pts <- size * pts %*% rot + rep(1, res) %*% t(c(x,y))
		polygon(pts, col=col, border=border, lwd=lwd)
	}

	# defaults
	if (missing(clades)) {
		tips <- phy$tip.label
		clades <- rep_len(1, length(tips))
		lines <- FALSE
		anc_labels <- FALSE
	} else if (!is.null(dim(clades))) {
		clades <- clades[match(phy$tip.label, clades[,1]),]
		tips <- clades[,1]
		clades <- clades[,2]
	} else if (is.vector(clades)) {
		tips <- phy$tip.label
		if (length(tips) != length(clades)) {
			stop("length of 'clades'must be equal to the number of tree tips")
		}
	}
	L <- sort(unique(clades))
	K <- length(L)
	N <- length(tips)
	if (missing(color)) {
		if (length(unique(clades)) == 1) {
			color <- setNames(stem, L)
		} else {
			palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")		
			color <- setNames(palette[seq(K)], L)	
		}
	} else if (is.null(names(color)) & length(color) >= K) {
		color <- setNames(color[seq(K)], L)
	} else if (any(names(color) %in% L)) {
		color <- setNames(color[L], L)
	} else if (any(names(color) %in% tips)) {
		color <- setNames(color[tips], tips)
	}
	if (is.numeric(color)) {
		color <- setNames(grDevices::palette()[color], names(color))
	}

	# data preparation
	species <- split(tips, clades)
	if (isTRUE(order)) {
		phy <- ape::rotateConstr(phy, unlist(species[order(names(species))]))
	} else if (is.character(order)) {
		phy <- ape::rotateConstr(phy, unlist(species[order]))		
	}
	phy <- ape::read.tree(text=ape::write.tree(phy))
	off <- lapply(species, getOffspring, phy=phy)
	anc <- lapply(species, getAncestor, phy=phy)
	paths <- mapply(getPaths, anc, off, MoreArgs=list(phy=phy), SIMPLIFY=FALSE)	
	group <- lapply(paths, makeEdges)
	group <- lapply(group, function(x) match(x[,2], phy$edge[,2]))
	group <- data.frame(edge=unlist(group), otu=rep(names(group), sapply(group, length)), stringsAsFactors=FALSE)
	stems <- setdiff(seq(nrow(phy$edge)), group$edge)
	if (length(stems) > 0) {
		stems <- data.frame(edge=stems, otu="stem", stringsAsFactors=FALSE)
	} else {
		stems <- NULL
	}	
	group <- rbind(group, stems)
	group <- group[order(as.numeric(group$edge)),2]
	tipcolor <- group
	nodes <- sort(unique(c(phy$edge[group == "stem",1], phy$edge[group != "stem",2])))
	stems <- nodes[sapply(nodes, function(i) any(group[phy$edge[,1] == i] != group[phy$edge[,2] == i]))]
	while (length(stems) > 0) {
		group[(phy$edge[,1] %in% stems | phy$edge[,2] %in% stems)] <- "stem"
		nodes <- sort(unique(c(phy$edge[group == "stem",1], phy$edge[group != "stem",2])))
		stems <- nodes[sapply(nodes, function(i) any(group[phy$edge[,1] == i] != group[phy$edge[,2] == i]))]
	}
	group <- ifelse(phy$edge[,2] <= N, tipcolor, group)
	if (!missing(edge.color)) {
		edge.col <- edge.color
	} else {
		edge.col <- c(color, stem=stem)[group]
	}

	# plotting
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "tree_clades.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}

# margins	
	if (missing(mai)) {
		mai <- c(0.02, 0.02, 0.02, 0.02)
		if (isTRUE(scale)) {
			cex.axis <- name.cex
			onechar <- strheight(s=0, units="inches", cex=1) * cex.axis
			mai[1] <- 2 * onechar
		}
		if (isTRUE(lines)) {
			strsize <- ifelse(name.las %in% c(1,2), strwidth, strheight)
			strw <- max(strsize(s=clades, units="inches", cex=1, font=name.font)) * name.cex
			onechar <- strsize(s="0", units="inches", cex=1, font=name.font) * name.cex
			mai[4] <- line.pos * 0.2 + (1.0 * onechar + strw + 0.5 * onechar) * isTRUE(names)
		}
	} else {
		mai <- rep_len(mai, 4)
	}

# layout
	xy <- ape::plotPhyloCoor(phy, type=type, direction="rightwards", use.edge.length=TRUE)
	ext <- 1 + 0.01 * extension
	xlim <- scale(range(xy[,1]), scale=FALSE)
	xlim <- as.numeric(ext * xlim + attr(xlim,"scaled:center"))
	if (isTRUE(tip_labels) | is.character(tip_labels)) {
		if (isTRUE(tip_labels)) {
			tip_labels <- phy$tip.label
		}
		strw <- strwidth(s=tip_labels, units="inches", font=tip.lab.font) * tip.lab.cex
		strw <- diff(xlim) * (strw / width) / (diff(c(mai[2],width-mai[4]))/width)
		xlim[2] <- max(xy[seq(N),1] + strw + 0.5 * strw / nchar(tip_labels) + tip.lab.pos * diff(xlim)) 
	}
	ylim <- scale(range(xy[,2]), scale=FALSE)
	ylim <- ext * ylim + attr(ylim,"scaled:center")

# plot
	suppressWarnings(par(mai=mai, new=TRUE))
	plot(xy, type="n", bty="n", xlim=xlim, ylim=ylim, axes=FALSE, ann=FALSE, xaxs="i", yaxs="i")
	drawTree(phy, xy, lwd=edge.width, col=edge.col)

# clade lines	
	if (isTRUE(lines)) {
		part <- ape::prop.part(phy)
		nodes <- sort(unique(phy$edge[group == "stem",2]))
		full <- nodes[sapply(nodes, function(i) sum(group[phy$edge[,1] == i] != "stem") == 2)]
		ranges <- lapply(part[full - N], range)
		half <- nodes[sapply(nodes, function(i) sum(group[phy$edge[,1] == i] != "stem") == 1)]
		if (length(half) > 0) {
			halves <- sapply(half, function(i) which(phy$edge[,1] == i)[group[phy$edge[,1] == i] != "stem"])
			halves <- as.list(phy$edge[halves,2])
			halves <- lapply(halves, function(i) if (i <= N) i else range(part[[i]]))
			ranges <- c(halves, ranges)
		}
		ranges <- ranges[order(sapply(ranges, min))]
		over <- line.over * strheight(s="0", units="user", cex=1, font=name.font) * name.cex
		for (i in which(sapply(ranges, length) == 1)) {
			ranges[[i]] <- ranges[[i]] + c(-1, 1) * over
		}
		col_axis <- color[group[match(as.integer(sapply(ranges, min)),phy$edge[,2])]]
		if (is.null(name.col)) {
			col_names <- col_axis
		} else {
			col_names <- rep_len(name.col, length(col_axis))
		}
		sp <- names(col_axis)
		for (i in seq_along(ranges)) {
			axis(side=4, at=ranges[[i]], labels=FALSE, tick=TRUE, lwd=edge.width, lwd.ticks=0, col=col_axis[i], line=line.pos)
		}
		if (isTRUE(names)) {
			for (i in seq_along(ranges)) {
				mtext(text=sp[i], side=4, at=mean(ranges[[i]]), line=1.0*onechar/0.2, cex=name.cex, font=name.font, las=name.las, col=col_names[i])
			}
		}
	}

# divergence scale	
	if (isTRUE(scale)) {
		at <- axTicks(side=1, usr=range(xy[,1]))
		step <- diff(at[1:2])
		at <- mean(at) + c(-0.5, 0.5) * step
		onechar <- strheight(s="0", units="inches", cex=1) * cex.axis			
		axis(side=1, at=at, labels=FALSE, lwd=edge.width/2, cex.axis=cex.axis, tick=TRUE, lwd.ticks=0)
		mtext(text=as.character(step), side=1, at=mean(at), cex=cex.axis, line=0.5*onechar/0.2)
	}

# tip labels	
	if (is.character(tip_labels)) {
		strw <- strwidth(tip_labels, units="inches", cex=1, font=1) * tip.lab.cex
		strw <- (strw - par("fig")[1]) / diff(par("fig")[1:2])
		strw <- (strw - mai[2]) / (width - mai[4] - mai[2]) * diff(xlim) + xlim[1]		
		tipxy <- xy[seq(N),]
		tipxy[,1] <- tipxy[,1] + tip.lab.pos * diff(range(xy[,1])) + 0.5 * strw
		text(tipxy, label=tip_labels, cex=tip.lab.cex, font=tip.lab.font)
	}

# node labels	
	if (isTRUE(node_labels)) {
		node_labels <- phy$node.label
	} else if (length(names(node_labels)) > 0) {
		node_labels <- node_labels[match(seq(phy$Nnode), as.numeric(names(node_labels)))]
	}
	if (!is.null(root_label)) {
		node_labels[1] <- root_label
	}
	if (!is.null(node_labels)) {
		node.lab.asp <- node.lab.asp * (diff(ylim) / diff(xlim))
		node.lab.asp <- rep_len(node.lab.asp, phy$Nnode)
		node.lab.cex <- rep_len(node.lab.cex, phy$Nnode)
		node.lab.pos <- rep_len(node.lab.pos, phy$Nnode)
		node.lab.col <- rep_len(node.lab.col, phy$Nnode)
		node.lab.bg <- rep_len(node.lab.bg, phy$Nnode)
		node.lab.lwd <- rep_len(node.lab.lwd, phy$Nnode)
		ntext.lab.cex <- rep_len(text.lab.cex, phy$Nnode)
		ntext.lab.font <- rep_len(text.lab.font, phy$Nnode)
		ntext.lab.col <- rep_len(text.lab.col, phy$Nnode)
		#strw <- max(strwidth(node_labels, units="user", cex=1, font=ntext.lab.font) * ntext.lab.cex)
		strw <- max(strwidth(node_labels, units="inches", cex=1, font=ntext.lab.font) * ntext.lab.cex)
		strw <- (strw - par("fig")[1]) / diff(par("fig")[1:2])
		strw <- (strw - mai[2]) / (width - mai[4] - mai[2]) * diff(xlim) + xlim[1]		
		nodesize <- 1.5 * strw
		nodexy <- xy[-seq(N),]
		if (any(node.lab.pos < 1)) {
			edge_lengths <- nodexy[phy$edge[,2],1] - nodexy[phy$edge[,1],1]
			nodexy[,1] <- nodexy[,1] - (1 - node.lab.pos) * edge_lengths
		}
		nodes <- which(nchar(node_labels) > 0 & !is.na(node_labels)) + N
		if (isTRUE(anc_labels)) {
			nodes <- intersect(nodes, unique(phy$edge[group == "stem",2]))
		}
		nodes <- nodes - N
		if (!is.null(root_label)) {
			nodes <- c(1, nodes)
		}
		for (i in nodes) {
			drawEllipse(x=nodexy[i,1], y=nodexy[i,2], size=nodesize, w=1, h=node.lab.asp[i], r=0, res=200, col=node.lab.bg[i], border=node.lab.col[i], lwd=node.lab.lwd[i])
		}
		text(nodexy[nodes,], labels=node_labels[nodes], col=ntext.lab.col[nodes], font=ntext.lab.font[nodes], cex=ntext.lab.cex[nodes], adj=text.lab.pos)
	}
	
# edge labels	
	if (is.character(edge_labels)) {
		if (length(names(edge_labels)) > 0) {
			edge_labels <- edge_labels[match(seq(nrow(phy$edge)), as.numeric(names(edge_labels)))]
		} else if (length(edge_labels) != nrow(phy$edge)) {
			warning("'edge_labels' must be named or equal in length to the number of group")
		}
		edges <- which(nchar(edge_labels) > 0 & !is.na(edge_labels))
		edge.lab.asp <- edge.lab.asp * (diff(ylim) / diff(xlim))
		edge.lab.asp <- rep_len(edge.lab.asp, length(edges))
		edge.lab.cex <- rep_len(edge.lab.cex, length(edges))
		edge.lab.pos <- rep_len(edge.lab.pos, length(edges))
		edge.lab.col <- rep_len(edge.lab.col, length(edges))
		edge.lab.bg <- rep_len(edge.lab.bg, length(edges))
		edge.lab.lwd <- rep_len(edge.lab.lwd, length(edges))
		etext.lab.cex <- rep_len(text.lab.cex, length(edges))
		etext.lab.font <- rep_len(text.lab.font, length(edges))
		etext.lab.col <- rep_len(text.lab.col, length(edges))
		strw <- max(strwidth(edge_labels, units="inches", cex=1, font=etext.lab.font) * etext.lab.cex)
		strw <- (strw - par("fig")[1]) / diff(par("fig")[1:2])
		strw <- (strw - mai[2]) / (width - mai[4] - mai[2]) * diff(xlim) + xlim[1]		
		edgesize <- 1.5 * strw
		edgexy <- xy[phy$edge[edges,2],]
		if (any(edge.lab.pos < 1)) {
			edge_lengths <- edgexy[phy$edge[,2],1] - edgexy[phy$edge[,1],1]
			edgexy[,1] <- edgexy[,1] - (1 - edge.lab.pos) * edge_lengths
		}
		for (i in edges) {
			drawEllipse(x=edgexy[i,1], y=edgexy[i,2], size=edgesize, w=1, h=edge.lab.asp[i], r=0, res=200, col=edge.lab.bg[i], border=edge.lab.col[i], lwd=edge.lab.lwd[i])			
		}
	}

# tidying up
	par(new=FALSE)
	
# closing pdf device	
	if (isTRUE(device == "pdf")) {
		dev.off()
	}
}



