readIpyrad <- function(file) {
	loci <- readLines(file)
	ends <- grep("//", loci)
	loci <- strsplit(loci[-ends], " +")
	nloci <- length(ends)
	nind <- indperloc <- diff(c(0,ends)) - 1
	indnames <- sapply(loci, "[", 1)
	indnames <- split(indnames, rep(seq(nloci), nind))
	loci <- sapply(loci, "[", 2)
	loci <- split(loci, rep(seq(nloci), nind))
	loci <- lapply(loci, function(x) do.call(rbind, strsplit(x, "")))
	for (i in seq(nloci)) {
		rownames(loci[[i]]) <- indnames[[i]]
	}
	return(loci)
}
