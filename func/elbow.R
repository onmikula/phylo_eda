# elbow method for finding of a breakpoint
# Salvador, S. and Chan, P. (2004) Determining the number of clusters/segments in hierarchical clustering/segmentation algorithms. In Proceedings of the Sixteenth IEEE International Conference on Tools with Artificial Intelligence, pp. 576–584, Institute of Electrical and Electronics Engineers. DOI: 10.1109/ICTAI.2004.50

# ARGUMENTS
# y - sekvence of evaluation metric (objective function) 
# x - rank (default) or spacing of evaluated solutions

elbow <- function(y, x=seq_along(y)) {
	k <- length(y)
	b <- 2:(k - 2)
	w <- rbind(b, (k-b)) / k
	lms <- lapply(b, function(p) list(L=lm(y[1:p] ~ x[1:p]), R=lm(y[(p+1):k] ~ x[(p+1):k])))
	rmse <- colSums(w * sqrt(sapply(lms, function(m) sapply(m, deviance))))
	result <- list(x=x, y=y, b=b, rmse=rmse, best=b[which.min(rmse)])
	class(result) <- "elbow"
	return(result)
}


plot.elbow <- function(elb) {
	par(mai=c(1.02, 1.02, 0.42, 0.42))
	plot(elb$x, elb$y, ylim=range(c(elb$y, elb$rmse)), xlab="breakpoints", ylab="RMSE", type="b", col="darkgrey", cex.lab=1.5, cex.axis=1.15, pch=21, cex=1.15, lwd=1.5)
	lines(elb$x[2:(length(elb$x)-2)], elb$rmse, type="b", pch=21, col="black", bg=ifelse(seq_along(elb$y) == (elb$best-1), "black", "white"), cex=1.15, lwd=1.5)
}

