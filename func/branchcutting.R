# depends on: ape

# inputs:
#	phy - object of class "phylo", treated as unrooted
#	method - method of breakpoint estimation
#	het, con - two-column matrices where each row contain labels of tips that are heterospecific or conspecific

branchcutting <- function(phy, method=c("step","elbow"), het, con, kernel="rectangular") {
	upp <- upper.tri(diag(length(phy$tip.label)))
	loss <- numeric(length(phy$edge.length))
	for (i in seq_along(phy$edge.length)){
		itre <- phy
		itre$edge.length[i] <- 0
		loss[i] <- mean(cophenetic(itre)[upp]) 
	}
	loss <- mean(cophenetic(phy)[upp]) - loss
	ord <- order(loss)
	loss <- loss[ord]

	method <- method[1]
	if (method == "step") {
		dif <- diff(loss)
		bw <- mean(as.numeric(dist(dif)))
		dens <- density(dif, bw=bw, adjust=1, kernel=kernel, cut=0)
		d <- sort(ifelse(dens$y < 0, 0, dens$y), decreasing=TRUE)
		mindens <- d[min(which(cumsum(d / sum(d)) >= 0.9999))]
		if (all(dens$y <= mindens)) {
			bpoint <- length(loss)
		} else {
			separpeaks <- which(diff(as.numeric(dens$y > mindens)) == 1)
			if (length(separpeaks) > 0) {
				thr <- dens$x[min(separpeaks) + 1]
				bpoint <- length(loss) - min(which(dif >= thr))
			} else {
				bpoint <- 0
			}
		}
	}
	if (method == "elbow") {
		tloss <- loss * (length(loss) - 1) / diff(range(loss))
		AB <- sqrt(diff(range(tloss))^2 + (length(tloss)-1)^2)
		AC <- sapply(seq_along(tloss), function(i) sqrt((i-1)^2 + (tloss[i]-min(tloss))^2))
		BC <- sapply(seq_along(tloss), function(i) sqrt((length(tloss) - i)^2 + (max(tloss)-tloss[i])^2))
		s = (AB + AC + BC) / 2
		a = sqrt(s * (s - AB) * (s - AC) * (s - BC))
		h = a / (0.5 * AB)
		bpoint <- length(tloss) - which.max(h)
	}

	if (!missing(het)) {
		if (is.vector(het)) {
				het <- matrix(het, 1, 2)
		}
		het <- matrix(match(het, phy$tip.label), nrow(het), 2)
		if (is.rooted(phy)) {
			nodes <- unlist(lapply(seq(nrow(het)), function(i) nodepath(phy, het[i,1], het[i,2])))
			hpoint <- match(which(phy$edge[,1] %in% nodes), ord)
		} else {
			nodes <- lapply(seq(nrow(het)), function(i) nodepath(phy, het[i,1], het[i,2]))
			branches <- lapply(nodes, function(x) cbind(x[-length(x)],x[-1]))
			branches <- do.call(rbind, lapply(branches, function(x) rbind(x,x[,2:1])))
			branches <- split(branches, seq(nrow(branches)))
			branches <- Map(rbind, branches, rep(list(phy$edge), length(branches)))
			branches <- unname(unlist(lapply(branches, function(x) which(duplicated(x)) - 1)))
			hpoint <- match(branches, ord)
		}
		hpoint <- length(loss) - hpoint[which.max(loss[hpoint])] + 1
	} else {
		hpoint <- 0
	}
	if (!missing(con)) {
		if (is.vector(con)) {
				con <- matrix(con, 1, 2)
		}
		con <- matrix(match(con, phy$tip.label), nrow(con), 2)
		if (is.rooted(phy)) {
			nodes <- unlist(lapply(seq(nrow(con)), function(i) nodepath(phy, con[i,1], con[i,2])))
			cpoint <- match(which(phy$edge[,1] %in% nodes), ord)
		} else {
			nodes <- lapply(seq(nrow(con)), function(i) nodepath(phy, con[i,1], con[i,2]))
			branches <- lapply(nodes, function(x) cbind(x[-length(x)],x[-1]))
			branches <- do.call(rbind, lapply(branches, function(x) rbind(x,x[,2:1])))
			branches <- split(branches, seq(nrow(branches)))
			branches <- Map(rbind, branches, rep(list(phy$edge), length(branches)))
			branches <- unname(unlist(lapply(branches, function(x) which(duplicated(x)) - 1)))
			cpoint <- match(branches, ord)
		}
		cpoint <- length(loss) - cpoint[which.max(loss[cpoint])]
	} else {
		cpoint <- length(loss)
	}
	if (hpoint >= cpoint) {
		warning("heterospecific and conspecific constraints are in conflict and could not be applied")
	} else {
		bpoint <- min(c(max(c(bpoint, hpoint)), cpoint))
	}

	if (bpoint > 0) {
		collapsed <- phy
		collapsed$edge.length <- ifelse(seq_along(ord) %in% tail(ord, bpoint), collapsed$edge.length, 0)
		retnodes <- unique(sort(collapsed$edge[collapsed$edge.length > 0,]))
		if (is.rooted(phy)) {
			retnodes <- c(retnodes, length(phy$tip.label) + 1)
		}
		retpaths <- vector("list", 0.5 * length(retnodes) * (length(retnodes)-1))
		ij <- 1
		for (i in 1:(length(retnodes)-1)) {
			for (j in (i+1):length(retnodes)) {
				retpaths[[ij]] <- nodepath(phy, retnodes[i], retnodes[j])
				ij <- ij + 1
			}
		}
		retnodes <- unique(unlist(retpaths))
		retedges <- collapsed$edge[,1] %in% retnodes & collapsed$edge[,2] %in% retnodes
		collapsed$edge.length[retedges] <- phy$edge.length[retedges]
		if (!missing(con)) {
			dst=cophenetic(collapsed)
			conflict <- which(apply(con, 1, function(x, dst) dst[x[1], x[2]] > 0, dst=cophenetic(collapsed)))
			if (length(conflict) > 0) {
				if (!is.rooted(phy)) {
					top <- phy
					top$edge.length <- rep(1, length(top$edge.length))
					top <- cophenetic(top)
					tips <- con[conflict,,drop=FALSE]
					add <- matrix(,0,2) 
					for (i in 1:nrow(tips)) {
						temp <- top[tips[i,],]
						temp <- temp[,-as.numeric(apply(dst[tips[i,],] == 0, 1, which))]
						temp <- colnames(temp)[which(temp == min(temp), arr.ind=TRUE)[,2]]
						add <- rbind(add, cbind(tips[i,], which(colnames(dst) == temp)))
					}
					conflict <- c(conflict, nrow(con) + 1:nrow(add))
					con <- rbind(con, add)
				}
				nodes <- lapply(conflict, function(i) nodepath(collapsed, con[i,1], con[i,2]))
				branches <- lapply(nodes, function(x) cbind(x[-length(x)],x[-1]))
				branches <- do.call(rbind, lapply(branches, function(x) rbind(x,x[,2:1])))
				branches <- split(branches, seq(nrow(branches)))
				branches <- Map(rbind, branches, rep(list(phy$edge), length(branches)))
				branches <- unname(unlist(lapply(branches, function(x) which(duplicated(x)) - 1)))
				coledges <- unique(branches[collapsed$edge.length[branches] > 0])
				if (is.rooted(phy)) {
					old <- 0
					new <- sum(collapsed$edge.length[coledges])
					while (new > old) {
						old <- new
						coledges <- unique(c(coledges, which(collapsed$edge[,1] %in% collapsed$edge[coledges,2])))
						new <- sum(collapsed$edge.length[coledges])
					}
				}
				collapsed$edge.length[coledges] <- 0
			}
		}
		otus <- unique(lapply(seq_along(phy$tip.label), function(i, dst) names(which(dst[i,] == 0)), dst=cophenetic(collapsed)))
	} else {
		collapsed <- phy
		otus <- list(phy$tip.label)
	}

	result <- list(phy=phy, branch=ord, loss=loss, bpoint=bpoint, collapsed=collapsed, otus=otus)
	attr(result, "method") <- method
	attr(result, "kernel") <- kernel
	attr(result, "constraints") <- c("heterospecific"[!missing(het)], "conspecific"[!missing(con)])
	if (length(attr(result, "constraints")) == 0 | hpoint >= cpoint)
		attr(result, "constraints") <- NULL
	class(result) <- "bcut"

	return(result)
}


plot.bcut <- function(bcut, which=c("tree","otus","loss","dens"), edge.color, edge.width, pt.cex, pt.col=1, ...) {
	if (missing("edge.color")) {
		edge.color <- 1
	}
	if (missing("edge.width")) {
		edge.width <- 1
	}
	if (missing("pt.cex")) {
		pt.cex <- 1
	}
	which <- which[1]
	if (which == "tree") {
		col <- (seq_along(bcut$branch) %in% rev(bcut$branch)[seq(bcut$bpoint)]) + 1
		col[col == 1 & bcut$collapsed$edge.length > 0] <- 3
		col[col == 2 & bcut$collapsed$edge.length == 0] <- 4
		if (is.rooted(bcut$phy)) {
			plot(bcut$phy, type="phylogram", edge.color=col, edge.width=edge.width, ...)
		} else {
			plot(bcut$phy, type="unrooted", edge.color=col, edge.width=edge.width, ...)
		}
	}
	if (which == "otus") {
		size <- sapply(bcut$otus, length)
		prop <- size / sum(size)
		retained <- which(bcut$collapsed$edge.length > 0)
		anc <- lapply(bcut$otus, function(x) getMRCA(bcut$phy, x))
		nz <- sort(unique(as.numeric(bcut$phy$edge[retained,])))
		dst <- dist.nodes(bcut$collapsed)
		for (i in seq_along(anc)) {
			if (size[i] == 1) {
				anc[[i]] <- match(bcut$otus[[i]], bcut$collapsed$tip.label)
			}
			if (!anc[[i]] %in% nz) {
				anc[[i]] <- nz[which(dst[anc[[i]],nz] == 0)]
			}
		}
		anc <- unlist(anc)
		type <- c("unrooted","cladogram")[is.rooted(bcut$phy) + 1]
		plot(bcut$phy, type=type, show.tip.label=FALSE, edge.color="transparent", edge.width=edge.width, ...)
		lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		xy <- cbind(lastPP$xx, lastPP$yy)
		xspl <- vector("list", length(anc))
		ultrametr <- is.ultrametric(bcut$phy, tol=0.0001)
		for (i in seq_along(anc)) {
			tips <- match(bcut$otus[[i]], bcut$phy$tip.label)
			nodes <- lapply(tips, function(x) nodepath(bcut$phy, x, anc[i]))
			nodes <- unique(unlist(nodes))
			nodes <- nodes[chull(xy[nodes,])]
			if (length(nodes) >= 3) {
				xynodes <- scale(xy[nodes,], scale=FALSE)
				xynodes <- 1.1 * xynodes + rep(1, length(nodes)) %*% t(attr(xynodes,"scaled:center"))
				spl <- c(-0.5, 0)[ultrametr + 1]
				xspl[[i]] <- xspline(x=xynodes[,1], y=xynodes[,2], s=rep(spl, length(nodes)), open=FALSE, draw=FALSE)
			}
		}
		if (ultrametr) {
			xmax <- max(unlist(lapply(xspl, "[[", "x")))
			xtip <- max(xy[Ntip(bcut$phy),1])
			for (i in which(!sapply(xspl, is.null))) {
				xspl[[i]]$x[xspl[[i]]$x >= xtip] <- xmax
			}
		}
		for (i in which(!sapply(xspl, is.null))) {
			x <- data.frame(x=xspl[[i]]$x, y=xspl[[i]]$y)
			ch <- grDevices::chull(x)
			xspl[[i]] <- x[c(ch,ch[1]),]
		}		
		for (i in which(!sapply(xspl, is.null))) {
			polygon(xspl[[i]], col="grey80", border="transparent")
		}
		for (i in seq_along(bcut$phy$edge.length)) {
			lines(xy[bcut$collapsed$edge[i,],], lwd=edge.width, col=edge.color)
		}
		points(xy[anc,,drop=FALSE], pch=16, cex=pt.cex, col=pt.col)
	}
	if (which == "loss") {
		plot(seq_along(bcut$loss), bcut$loss, type="n", xlab="rank", ylab="loss", ...)
		if (is.null(bcut$bpoint))
			ind <- seq_along(bcut$loss)
		else
			ind <- seq_along(bcut$loss) < (length(bcut$loss) - bcut$bpoint + 1)
		col <- ifelse(ind, "black", "red")
		col[col == "black" & bcut$collapsed$edge.length[bcut$branch] > 0] <- "green3"
		col[col == "red" & bcut$collapsed$edge.length[bcut$branch] == 0] <- "blue"
		subs <- col == "green3"
		points(seq_along(bcut$loss)[!subs], bcut$loss[!subs], pch=16, col=col[!subs], ...)
		points(seq_along(bcut$loss)[subs], bcut$loss[subs], pch=16, col=col[subs], ...)
	}
	if (which == "dens") {
		dif <- diff(bcut$loss)
		bw <- mean(as.numeric(dist(dif)))
		dens <- density(dif, bw=bw, adjust=1, kernel=attr(bcut, "kernel"), cut=0)
		d <- sort(ifelse(dens$y < 0, 0, dens$y), decreasing=TRUE)
		mindens <- d[min(which(cumsum(d / sum(d)) >= 0.9999))]
		thr <- dens$x[min(which(diff(as.numeric(dens$y > mindens)) == 1)) + 1]
		plot(dens, ...)
		abline(v=thr, col="red")
	}
}


summary.bcut <- function(bcut) {
	constraints <- attr(bcut, "constraints")
	if (is.null(constraints))
		constraints <- "none imposed"
	else if (identical(constraints, c("heterospecific", "conspecific")))
		constraints <- "heterospecific & conspecific"
	cat("Object of class 'bcut':\n")
	cat(paste(paste("\tNo. of OTUs:", length(bcut$otus)), "\n", sep=""))
	cat(paste(paste("\tNo. of tree tips:", length(bcut[["phy"]]$tip.label)), "\n", sep=""))
	cat(paste(paste("\tNo. of retained branches:", bcut[["bpoint"]]), "\n", sep=""))
	cat(paste(paste("\tMethod:", attr(bcut, "method")), "\n", sep=""))
	cat(paste(paste("\tConstraints:", constraints), "\n", sep=""))
}

							
print.bcut <- function(bcut) {
	k <- length(bcut$otus)
	nn <- sapply(bcut$otus, length)
	o <- ifelse(k == 1, "OTU", "OTUs")

	cat(paste(paste("Object of class 'bcut' with", k, "delimited", o), ":\n", sep=""))
	for (i in seq(k)) {
		six <- head(bcut$otus[[i]])
		if (nn[i] > 6)
			six <- c(six, "...")
		otu <- paste(paste("\t", paste("OTU", i), ":", sep=""), paste(six, collapse=", "))
		cat(paste(otu, "\n", sep=""))
	}
}


write.bcut <- function(bcut, file=NULL, return=TRUE) {
	df <- data.frame(ID=unlist(bcut$otus), OTU=rep(seq_along(bcut$otus), sapply(bcut$otus, length)), stringsAsFactors=FALSE)
	if (!is.null(file)) {
		write.table(df, file, row.names=FALSE, quote=FALSE, sep="\t")
	}
	if (isTRUE(return)) {
		attr(df, "nsp") <- length(bcut$otus)
		return(df)
	}
}

