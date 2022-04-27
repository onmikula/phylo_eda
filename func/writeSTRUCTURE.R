read_str_data <- function(x, row.names=1) {
	str <- read.delim(x, skip=1, header=FALSE, row.names=row.names)
	colnams <- paste(rep(unlist(strsplit(readLines(x, n=1), "[[:blank:]]+")), each=2), c(0,1), sep="_")
	add <- ncol(str) - length(colnams)
	if (add > 0) {
		atr <- str[,seq(add)]
		str <- str[,-seq(add)]
	} else {
		atr <- NULL
	}
	colnames(str) <- colnams
	#return(str[,grepl("^[-]*[[:digit:]]", str[1,])])
	attr(str, "info") <- atr
	return(str)
}


# it supposes phased diploid genotypes labelled by individual name + allele identifier specified by the argument 'allele'

write_str_data <- function(loci, file=NULL, return=FALSE, allele="_[[:digit:]]+$", linkage=FALSE, PopData=NULL, PopFlag=1, structure=NULL) {
if (is.null(structure)) {
	if (inherits(loci, "matrix")) {
		locnam <- colnames(loci)
		if (is.null(locnam)) {
			locnam <- paste0("L", formatC(seq(ncol(loci)), format="d", flag=0, width=nchar(ncol(loci))))
		}
		loci <- setNames(lapply(seq(ncol(loci)), function(i) loci[,i,drop=FALSE]), locnam)
	} else {
		if (is.null(names(loci))) {
			names(loci) <- paste0("L", formatC(seq_along(loci), format="d", flag=0, width=nchar(length(loci))))
		}
		for (i in seq_along(loci)) {
			if (is.null(colnames(loci[[i]])) & ncol(loci[[i]]) > 1) {
				colnames(loci[[i]]) <- paste0(names(loci)[i], "S", formatC(seq(ncol(loci[[i]])), format="d", flag=0, width=nchar(ncol(loci[[i]]))))
			} else if (is.null(colnames(loci[[i]])) & ncol(loci[[i]]) == 1) {
				colnames(loci[[i]]) <- names(loci)[i]
			}
		}
	}
	char <- Filter(function(x) sign(x) == 1, suppressWarnings(as.numeric(loci[[1]])))
	type <- ifelse(length(char) == 0, "snps", "num")
	for (i in seq_along(loci)) {
		loci[[i]] <- loci[[i]][order(rownames(loci[[i]])),,drop=FALSE]
		if (type == "snps") {
			loci[[i]] <- toupper(loci[[i]])
			#loci[[i]][loci[[i]] %in% c("A","C","G","T")] <- NA
			for (j in seq(ncol(loci[[i]]))) {
				loci[[i]][,j] <- match(loci[[i]][,j], c("A","C","G","T"))
			}
		}
		ind <- gsub(allele, "", rownames(loci[[i]]))
		loc <- paste(rep(colnames(loci[[i]]), each=2), c(0,1), sep="_")
		loci[[i]] <- matrix(unlist(split(loci[[i]], ind)), length(ind)/2, length(loc), byrow=TRUE, dimnames=list(unique(ind), loc))
	}
	rows <- sort(unique(unlist(lapply(loci, rownames))))
	structure <- vector("list", length(loci))
	for (i in seq_along(structure)) {
		structure[[i]] <- matrix(-9, length(rows), ncol(loci[[i]]))
		structure[[i]] <- loci[[i]][match(rows, rownames(loci[[i]])),]
		structure[[i]][,] <- match(structure[[i]], sort(na.omit(unique(as.vector(structure[[i]])))))
		dimnames(structure[[i]]) <- list(rows, colnames(loci[[i]]))
	}
	if (isTRUE(linkage)) {
		linkage <- paste(unlist(lapply(structure, function(x) rep(c(-1,0), c(1, ncol(x)/2 - 1)))),  collapse="\t")
	}
	structure <- do.call(cbind, structure)
	structure[is.na(structure)] <- "-9"
}
	if (isFALSE(linkage)) {
		linkage <- NULL
	} else {
		linkage <- paste(linkage,  collapse="\t")
	}
	if (!is.null(PopData)) {
		matching <- match(rownames(structure), PopData[,1])
		PopData <- PopData[matching,2]
		PopFlag <- rep_len(as.numeric(PopFlag), length(PopData))[matching]
	} else {
		PopFlag <- NULL
	}
	locnames <- paste(unique(gsub("_[[:digit:]]+$", "", colnames(structure))),  collapse="\t")
	header <- cbind(rownames(structure), PopData, PopFlag)
	input <- cbind(header, structure)
	input <- apply(input, 1, paste,  collapse="\t")
	input <- c(locnames, linkage, input)
	if (!is.null(file)) {
		writeLines(input, con=file)
	}
	if (isTRUE(return)) {
		return(structure)
	}
}






parse_str_params <- function(df, templates =c("mainparams", "extraparams"), outputs=c("main.txt", "extra.txt")) {
	names(df) <- c("Arg", "Value")
	df$Arg <- as.character(paste0("#define[[:blank:]]+", toupper(df$Arg)))

	main <- readLines(templates[1])
	hits <- lapply(df$Arg, grep, x=main)
	mainpar <- which(sapply(hits, length) == 1)
	for (i in mainpar) {
		legend <- sub("^#[[:alnum:][:blank:]]+\\/", "/", main[hits[[i]]])
		main[hits[[i]]] <- paste(sub("[[:blank:]]+", " ", df$Arg[i], fixed=TRUE), "   ", df$Value[i], "   ", legend)
	}
	writeLines(main, con=outputs[1])
	df <- df[-mainpar,,drop=FALSE]
	extra <- readLines(templates[2])
	if (nrow(df) > 0) {
		hits <- lapply(df$Arg, grep, x=extra)
		extrapar <- which(sapply(hits, length) == 1)
		for (i in extrapar) {
			legend <- sub("^#[[:alnum:][:blank:]]+\\/", "/", extra[hits[[i]]])
			extra[hits[[i]]] <- paste(sub("[[:blank:]]+", " ", df$Arg[i], fixed=TRUE), "   ", df$Value[i], "   ", legend)
		}
	}
	writeLines(extra, con=outputs[2])	
}


### read_str_output ###
# file: name of the file with output of STRUCTURE (so far it works for admixture model with no POPINFO prior)

read_str_output <- function(file) {
	x <- suppressWarnings(readLines(file))
	x <- x[nchar(x) > 0]
	K <- as.numeric(sub("MAXPOPS=", "", regmatches(x, regexpr("MAXPOPS=[[:digit:]]+", x))))
	nind <- as.numeric(sub("NUMINDS=", "", regmatches(x, regexpr("NUMINDS=[[:digit:]]+", x))))
	start <- grep("Inferred ancestry of individuals", x) + 2
	end <- start + nind - 1
	qmat <- x[start:end]
	qmat <- gsub("^[[:blank:]]+", "", qmat)
	qmat <- gsub(":", " ", qmat)	
	qmat <- do.call(rbind, strsplit(qmat, "[[:blank:]]+"))
	qdf <- data.frame(Label=qmat[,2], Miss=as.numeric(gsub("[\\(\\)]", "", qmat[,3])))
	est <- qmat[,3+(1:K)]
	mode(est) <- "numeric"
	rownames(est) <- qdf$Label
	colnames(est) <- paste0("Cluster", formatC(seq(K), format="d", flag=0, width=nchar(K)))
	for (i in seq(K)) {
		int <- gsub("[\\(\\)]", "", qmat[,K + 3 + i])
		int <- do.call(rbind, lapply(strsplit(int, ","), as.numeric))
		colnames(int) <- paste(colnames(est)[i], c("low","upp"), sep="_")
		est <- cbind(est, int)
	}
	est <- as.data.frame(est)
	attr(est, "K") <- K
	attr(est, "Missing") <- qdf$Miss
	return(est)
}


### plot_str_output ###
# out: data frame produced by read_str_output
# palette: colors assigned to the clusters
# ord: order of individuals in the plot, if TRUE (default) it is determined internally, if FALSE it is ignored 
# device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
# file: name of pdf file if device == "pdf"
# width: width of the device in inches
# height: height of the device in inches
# mai: 
# interval: whether to display width of probability interval by color transparency
# expr: expression allowing to add other elements (e.g. legend)

plot_str_output <- function(out, palette, ord=TRUE, device, file, width=10, height=5, mai=NULL, interval=FALSE, expr) {
	N <- nrow(out)
	K <- max(as.numeric(gsub("Cluster", "", unlist(regmatches(names(out), gregexpr("Cluster[[:digit:]]+", names(out)))))))
	est <- out[,paste0("Cluster", formatC(seq(K), format="d", flag=0, width=nchar(K)))]
	low <- out[,paste0("Cluster", formatC(seq(K), format="d", flag=0, width=nchar(K)), "_low")]
	upp <- out[,paste0("Cluster", formatC(seq(K), format="d", flag=0, width=nchar(K)), "_upp")]
	alpha <- setNames(1 - (upp - low), names(est))
	if (!identical(sort(rownames(out)), sort(ord))) {
		if (isTRUE(ord)) {
			ord <- hclust(dist(est), method="average")$order
		} else if (isFALSE(ord)) {
			ord <- seq(N)
		} else {
			ord <- seq(N)
			warning("'ord' argument was set to FALSE")
		}	
	}	
	nam <- rownames(out)[ord]
	est <- est[ord,]
	alpha <- alpha[ord,]
	csum <- cbind(0, t(apply(est, 1, cumsum)))
	csum <- diag(1/csum[,K+1]) %*% csum

	if (missing(palette)) {
		palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080", White="#FFFFFF", Black="#000000")
	}
	color <- setNames(palette[seq(K)], names(est))
	if (isFALSE(interval)) {
		alpha[,] <- 1
	}
	color <- do.call(cbind, lapply(seq(K), function(i) opacity(color[i], alpha[,i])))
	if (missing(mai)) {
		mai <- c(0.42, 0.82, 0.42, 0.42)
	}

	if (missing(device)) {
		if (.Platform$OS.type == "unix") {
			device <- "quartz"
		} else {
			device <- "x11"
		}
	}
	if (isTRUE(device == "pdf")) {
		file <- ifelse(missing(file), "str_barplot.pdf", file)
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}
	par(mai=mai)
	plot(0, 0, xlim=c(0, nrow(out)), ylim=c(0,1), xaxs="i", yaxs="i", xaxt="n", xlab="", ylab="Ancestry proportions", type="n")
	for (i in seq(N)) {
		x <- rep(i, 5) - c(1,0,0,1,1)
		for (j in seq(K)) {
			xy <- cbind(x, csum[i,j+c(0,0,1,1,0)])
			polygon(xy, border=color[i,j], col=color[i,j])
		}
	}
	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		}
	}
	if (isTRUE(device == "pdf")) {
		dev.off()
	}
}



### opacity ###
# based on add_trans at: https://rdrr.io/github/mcsiple/mmrefpoints/man/add_trans.html
# Arguments:
# color - input color
# alpha - degree of opacity (from 0 to 1)	
# Details:
# this function defines opacity of a color, 0 being fully transparent and 1 being fully visable
# works with either color and opacity a vector of equal length, or one of the two of length 1.

opacity <- function(color, alpha) {
	num2hex <- function(x) {
		hex <- unlist(strsplit("0123456789ABCDEF",split=""))
		return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
	}
	if (length(color) != length(alpha) & !any(c(length(color), length(alpha)) == 1)) stop("Vector lengths not correct")
	if (length(color) == 1 & length(alpha) > 1) color <- rep(color, length(alpha))
	if (length(alpha) == 1 & length(color) > 1) alpha <- rep(alpha, length(color))
	alpha <- as.integer(255 * alpha)
	rgb <- rbind(grDevices::col2rgb(color), alpha)
	res <- paste("#", apply(apply(rgb, 2, num2hex), 2, paste, collapse=""), sep="")
	return(res)	
}



writeGeno <- function(loci, file=NULL, return=FALSE) {
    alleles <- vector("list", length(loci))
    for (i in seq_along(alleles)) {
        bp <- matrix(toupper(loci[[i]]) %in% c("A","C","G","T"), nrow(loci[[i]]), ncol(loci[[i]]))
        sq <- loci[[i]][,apply(bp, 2, all), drop=FALSE]
        sq <- apply(sq, 1, paste, collapse="")
        sq <- sq[order(rownames(loci[[i]]))]
        ind <- gsub("_[[:digit:]]+$", "", rownames(loci[[i]]))
        alleles[[i]] <- setNames(match(sq, unique(sq)), ind)
    }
    rows <- sort(unique(unlist(lapply(alleles, names))))
    geno <- vector("list", length(alleles))
    for (i in seq_along(alleles)) {
        tab <- table(alleles[[i]], names(alleles[[i]]))
        geno[[i]] <- matrix(9, nrow(tab), length(rows), dimnames=list(NULL, rows))
        geno[[i]][, match(colnames(tab), rows)] <- tab
    }
    geno <- do.call(rbind, geno)
    geno <- apply(geno, 1, paste, collapse="")
    if (!is.null(file)) {
        writeLines(geno, con=file)
    }
    if (isTRUE(return)) {
        return(geno)
    }
}

