readStacks <- function(file, n=-1L, gcoords=FALSE) {

	X <- readLines(file, n=n)[-1]

    nams <- grep("^>", X)
    loci <- as.numeric(sapply(strsplit(X[nams], "_"), "[", 2))
    loci <- paste("Locus", formatC(loci, format="d", flag=0, width=max(nchar(loci))), sep="_")
    L <- unique(loci)
    if (isTRUE(gcoords)) {
        G <- X[nams][match(L, loci)]
        G <- regmatches(G, gregexpr(";.+]$", G))
        G <- gsub("^;\\s|]$", "", unlist(G))
        G <- do.call(rbind, setNames(strsplit(G, ",\\s"), L))
        G <- setNames(as.data.frame(G[,c(1,3,2)]), c("Scaffold", "Thread", "First"))
        mode(G$First) <- "numeric"
    }
    X <- setNames(X[-nams], gsub("^\\[|;$", "", regmatches(X[nams], regexpr("\\[.+;", X[nams]))))
    names(X) <- paste(names(X), paste(rep("A", length(X)), c(0,1), sep=""), sep="_")
    X <- split(X, loci)
    for (i in seq_along(X)) {
        X[[i]] <- ape::as.DNAbin(do.call(rbind, strsplit(X[[i]], "")))
        #if (i %% 1000 == 0) print(i)
    }
    if (isTRUE(gcoords)) {
    		G <- G[match(names(X), rownames(G)),]
        G$Last <- G$First + sapply(X, ncol)
        attr(X, "gcoords") <- G
    }
    
    return(X)

}

