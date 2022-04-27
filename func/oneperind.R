# loci - list of alignments with phased alleles (two per individual)
# allele - pattern (regular expression) of allele specifier in sequence names (default = last underscore followed by one or more alphanumeric characters)
# remove <- whether to remov allele specifier from sequence names

oneperind <- function(loci, allele="_[[:digit:]]$", remove=FALSE) {
	alleles <- sort(unique(unlist(lapply(loci, rownames))))
	indiv <- gsub(allele, "", alleles)
	for (i in seq_along(loci)) {
		als <- apply(loci[[i]], 1, paste, collapse="")
		als <- match(als, unique(als))
		ind <- indiv[match(rownames(loci[[i]]), alleles)]
		als_ind <- split(als, ind)
		het <- which(sapply(lapply(als_ind, unique), length) == 2)
		als_het <- setdiff(seq(max(als)), unlist(als_ind[-het]))
		tab <- table(als[als %in% als_het], ind[als %in% als_het])
		while (nrow(tab) > 0) {
			grp <- apply(tab, 1, paste, collapse="-")
			grp <- setNames(match(grp, unique(grp)), names(grp))
			dupl <- table(grp)
			dupl <- as.numeric(names(dupl)[dupl > 1])
			if (length(dupl) == 0) {
				dupl <- seq(nrow(tab))
			}
			for (j in dupl) {
				jallele <- intersect(rownames(tab), names(grp[grp==j]))
				if (length(jallele) > 0) {
					jallele <- sample(jallele, 1)
					jindiv <- colnames(tab)[tab[jallele,] != 0]
					jindiv <- sample(jindiv, length(jindiv))
					jindiv <- jindiv[which.min(colSums(tab[,jindiv,drop=FALSE]))]
					als_ind[[jindiv]] <- rep(as.numeric(jallele), 2)
					tab <- tab[rownames(tab) != jallele,,drop=FALSE]
					tab <- tab[,colnames(tab) != jindiv,drop=FALSE]
					tab <- tab[rowSums(tab) > 0,,drop=FALSE]
				} else {
					break
				}
			}
		}
		ii <- lapply(sapply(als_ind, length), sample, size=1)
		als_ind <- unlist(Map("[", als_ind, ii))
		als_nam <- unlist(Map("[", split(rownames(loci[[i]]), ind), ii))
		loci[[i]] <- loci[[i]][match(seq(max(als)), als),,drop=FALSE][als_ind ,,drop=FALSE]		
		if (isTRUE(remove)) {
			rownames(loci[[i]]) <- names(als_ind)		
		} else {
			rownames(loci[[i]]) <- als_nam
		}
	}
	return(loci)
}
