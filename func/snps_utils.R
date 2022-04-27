find_snps <- function(x, nalleles=4) {
	searching <- function(x, nalleles) {
		if (length(x) == 0) {
			return(c(var=0, pis=0))
		} else {
			x <- toupper(as.matrix(x))
			x[!x %in% c("A", "C", "G", "T") ] <- NA
			ndist <- apply(x, 2, dplyr::n_distinct, na.rm=TRUE)
			vs <- which(ndist > 1 & ndist <= nalleles)
			is <- which(sapply(vs, function(i) sum(table(x[, i]) > 1)) > 1)
			return(list(var=vs, pis=vs[is]))
		}
	}
	if (is.list(x)) {
		snps <- lapply(x, searching, nalleles=nalleles)
	} else {
		snps <- searching(x, nalleles=nalleles)
	}
	return(snps)
}


count_snps <- function(snps) {
	if (length(snps) == 2 & identical(names(snps), c("var","pis"))) {
		snps <- list(snps)
	}
	return(lapply(snps, function(x) lapply(x, length)))
}


get_snps <- function(loci, nalleles=4, type="pis", unique=FALSE, as.matrix=TRUE) {
	snps_pos <- lapply(find_snps(loci, nalleles= nalleles), "[[", type)
	nz <- sapply(snps_pos, length) > 0
	loci <- loci[nz]
	snps_pos <- snps_pos[nz]
	if (isTRUE(unique)) {
		snps_pos <- lapply(snps_pos, function(x) ifelse(length(x) == 1, x, sample(x, 1)))
	}
	snps <- Map(function(x, cols) x[,cols, drop=FALSE], loci, snps_pos)
	snps <- lapply(lapply(lapply(snps, as.matrix), as.character), toupper)	
	if (as.matrix == TRUE) {
		rows <- sort(unique(unlist(lapply(snps, rownames))))
		snps <- do.call(cbind, lapply(snps, function(x, r) x[match(r, rownames(x)),], r=rows))
		snps[is.na(snps)] <- "N"
	}
	return(snps)
}


make_binary_snps <- function(x, type=c("snps","structure"), center=FALSE, onerowperind=TRUE, allele="_[[:digit:]]$") {
	rows <- rownames(x)
	if (type[1] == "structure") {
		x[x %in% c("N","-","-9")] <- NA
		loci <- rep(1:(ncol(x)/2), each=2)
		for (i in 1:(ncol(x)/2)) {
			xi <- as.vector(x[,loci == i])
			x[,loci == i] <- match(xi, sort(unique(xi))) - 1
		}
	} else if (type[1] == "snps") {
		x[!x %in% c("A", "C", "G", "T") ] <- NA
		x <- apply(x, 2, function(n) match(n, sort(unique(n)))) - 1
	}
	mode(x) <- "numeric"
	rownames(x) <- rows
	if (isTRUE(onerowperind) & all(grepl(allele, rownames(x)))) {
		x <- do.call(rbind, split(x, gsub(allele, "", rownames(x))))
	}
	if (isFALSE(onerowperind) & all(!grepl(allele, rownames(x)))) {}
	if (isTRUE(center)) {
		x <- 2 * (x - 0.5)
		x[is.na(x)] <- 0 
	}
	return(x)
}



count_binary_snps <- function(b) {
	center <- any(b == -1)
	if (center) {
		b <- b / 2 + 0.5
		b[b == 0.5] <- NA
	}
	nloci <- ncol(b) / 2
	b <- do.call(cbind, by(t(b), rep(seq(nloci), each=2), colSums))
	if (center) {
		b <- 2 * (b - 1)
		b[is.na(b)] <- 0 
	}
	return(b)
}



bpcoverage <- function(x) {
	x[!toupper(x) %in% c("A", "C", "G", "T") ] <- NA
	return(colSums(!is.na(x)) / nrow(x))
}



indcoverage <- function(x) {
	seqnam <- sort(unique(unname(unlist(lapply(x, rownames)))))
	return(sapply(x, nrow) / length(seqnam))
}




