# reads files produced by mPTP
read_ptp <- function(file) {
	ptp <- readLines(file)
	df <- grep("^Species", ptp)
	nsp <- length(df)
	df <- cbind(df + 1, c(df[-1] - 2, length(ptp)))
	df <- lapply(seq(nsp), function(i, delim) delim[df[i,1]:df[i,2]], delim=ptp)
	df <- data.frame(ID=unlist(df), OTU=rep(seq(nsp), sapply(df, length)), stringsAsFactors=FALSE)
	attr(df, "nsp") <- nsp
	return(df)
}


