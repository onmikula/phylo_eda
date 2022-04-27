match_otus <- function(R, S, col.names=NULL) {
	sort_match <- function(x) setdiff(order(x, decreasing=TRUE), which(x == 0))
	ind <- sort(intersect(R[,1], S[,1]))
	if (any(c(nrow(R), nrow(S)) != length(ind))) {
		warning("matching limited to the common subset of individuals")
	}
	R <- R[match(ind, R[,1]),2]
	S <- S[match(ind, S[,1]),2]
	tab <- table(R, S)
	ord <- unique(unlist(apply(tab, 1, sort_match)))
	otus <- data.frame(ID=ind, OTU1=R, OTU2=ord[S], stringsAsFactors=FALSE)
	if (!is.null(col.names)) {
		if (length(col.names) == 2) {
			names(otus)[2:3] <- col.names
		} else {
			names(otus) <- col.names[1:3]
		}
	}
	return(otus)
}
