concatenate <- function(dna, as.DNAbin=FALSE, missing="-") {
	bp <- cumsum(sapply(dna, ncol))
	part <- attr(dna, "part")
	if (is.null(part)) {
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rownames(part) <- names(bp)
	}
	rows <- sort(unique(unlist(lapply(dna, rownames))))
	concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
	for (i in seq_along(dna)) {
		concatenated[rownames(dna[[i]]), part[i,1]:part[i,2]] <- dna[[i]]
	}
	if (isTRUE(as.DNAbin)) {
		concatenated <- ape::as.DNAbin(concatenated)
	}
	attr(concatenated, "part") <- part
	return(concatenated)
}
