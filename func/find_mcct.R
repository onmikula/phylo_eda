### FUNCTION
# find_mcct
### DESCRIPTION
# employs 'maxCladeCred' of 'phangorn' (Schliep 2011) to find maximum clade credibility tree with specified type of node heights (i.e., branch lengths)
### ARGUMENTS
# trees: object of class 'multiPhylo' or list of 'phylo' objects
# file: name of file to write the tree in
# return: whether to return the tree as 'phylo' object
# method: method to calculate node heights (assumed to be 'mean' or 'median'), if 'NA' or 'topology' is specified, the tree has no branch lengths
# monophyletic: whether to calculate node height only from monophyletic clades 
# target: target tree (as 'phylo' object) upon which the node heights from 'trees' are projected (instead of on the MCC tree estimated from them)
# burnin: burn-in proportion of the posterior sample specified as a porportion, i.e. interval [0,1), or a percentage, interval [0, 100)
# thinning: amount of subsampling of 'trees', specified as a proportion (if < 1) or as every i-th tree to be retained (if > 1)
# digits: the precision of posterior probability estimates (the desired number of digits after the decimal point)
### VALUE
# object of class 'phylo' and/or file with the tree in newick format
### DETAILS
# largely emulates functions of TreeAnnotator of BEAST 2 (Bouckaert et al. 2014), but it does not assume the trees to be neither ultrametric nor rooted
# it returns the maximum clade credibility tree (Drummond & Bouckaert 2015, p. 162) or the target tree with a specified type of node heights estimated form a sample of trees (presumably posterior sample from a Bayesian analysis)
# the 'common ancestor heights' option (Drummond & Bouckaert 2015, p. 92) corresponds to 'monophyletic=FALSE', but it can be combined with both mean and median
# currently, it does not support calcalation of the highest posterior density intervals and 
# except for 'phangorn' (Schliep 2011) it depends on 'ape' (Paradis & Schliep 2019) 

find_mcct <- function(trees, file=NULL, return=TRUE, method="mean", monophyletic=FALSE, target=NULL, burnin=0, thinning=1, digits=2) {
	if (burnin > 0) {
		burnin <- ifelse(burnin >= 1, 0.01 * burnin, burnin)
		trees <- trees[-seq(ceiling(burnin * length(trees)))]
	}
	if (thinning < 1) {
		trees <- trees[round(seq(1, length(trees), length=thinning * length(trees)))]
	} else if (thinning > 1) {
		trees <- trees[seq(1, length(trees), by=thinning)]
	}
	class(trees) <- "multiPhylo"
	rooted <- ape::is.rooted(trees[[1]])
	if (is.null(target)) {
		mcct <- phangorn::maxCladeCred(trees, tree=TRUE, part=NULL, rooted=rooted)	
	} else {
		if (is.character(target)) {
			target <- ape::read.tree(text=target)
		}
		if (!identical(rooted, ape::is.rooted(rooted))) {
			stop("both trees in the posterior sample and the target tree have to either rooted or unrooted")
		}
		mcct <- target
	}
	if (method != "topology" & !is.na(method)) {
		if (isTRUE(rooted)) {
			# HPD
			mcct <- node_heights(phy=mcct, trees=trees, method=method, monophyletic=monophyletic)	
		} else {
			mcct <- branch_lengths(phy=mcct, trees=trees, method=method)	
		}
	} else {
		mcct$edge.length <- NULL
	}
	digits <- min(digits, nchar(length(trees)))
	pp <- prob_clades(phy=mcct, psample=trees)
	mcct$node.label <- formatC(pp, format="f", digits=digits)
	if (!is.null(file)) {
		ape::write.tree(mcct, file)
	}
	if (isTRUE(return)) {
		return(mcct)
	}
}


### FUNCTION
# node_heights
### DESCRIPTION
# estimates specified type of node heights in a tree from a sample of trees
### ARGUMENTS
# phy: tree (class 'phylo') whose node heights are to be estimated
# trees: sample (class 'multiPhylo' or list of 'phylo' objects) providing information about the node heights
# method: method to calculate node heights (assumed to be 'mean' or 'median'), if 'NA' or 'topology' is specified, the tree has no branch lengths
# monophyletic: whether to calculate node height only from monophyletic clades 
### VALUE
# object of class 'phylo' and/or file with the tree in newick format

node_heights <- function(phy, trees, method="mean", monophyletic=FALSE) {
	fun <- match.fun(method)
	R <- ape::Ntip(phy) + 1
	tips <- phy$tip.label
	torders <- lapply(lapply(trees, "[[", "tip.label"), order)
	theights <- matrix(, length(tips), length(trees))
	for (i in seq_along(trees)) {
		tpaths <- ape::nodepath(trees[[i]])
		tbranches <- lapply(tpaths, function(x) match(x, trees[[i]]$edge[,2]))
		theights[,i] <- sapply(tbranches, function(x) sum(trees[[i]]$edge.length[x], na.rm=TRUE))[torders[[i]]]
	}
	theights <- theights[match(tips, sort(tips)),]
	TH <- apply(theights, 1, fun)
#	TH <- TH[match(tips, sort(tips))]
	clades <- lapply(ape::prop.part(phy), function(x) tips[x])
	nheights <- matrix(, length(clades), length(trees))
	for (i in seq_along(trees)) {
		if (isTRUE(monophyletic)) {
			monophyl <- which(sapply(clades, function(x) ape::is.monophyletic(trees[[i]], x)))
		} else {
			monophyl <- rep_len(TRUE, length(clades))
		}
		anc <- sapply(clades[monophyl], function(x) ape::getMRCA(trees[[i]], x))
		npaths <- lapply(anc, function(a) ape::nodepath(trees[[i]], from=R, to=a))
		nbranches <- lapply(npaths, function(x) match(x, trees[[i]]$edge[,2]))
		nheights[monophyl,i] <- sapply(nbranches, function(x) sum(trees[[i]]$edge.length[x], na.rm=TRUE))
	}
	NH <- apply(nheights, 1, fun, na.rm=TRUE)
	H <- c(TH, NH)
	phy$edge.length <- as.numeric(diff(t(matrix(H[phy$edge], ape::Nedge(phy), 2))))
	return(phy)
}


### FUNCTION
#branch_lengths
branch_lengths <- function(phy, trees, method="mean", monophyletic=FALSE) {	
	splits <- function(tre, tips) {
		n <- ape::Ntip(tre)
		spl <- vector("list", ape::Nedge(tre))
		spl[match(seq(n), tre$edge[,2])] <- as.list(tre$tip.label)
		spl[tre$edge[,2] > n] <- lapply(setdiff(tre$edge[,2], seq(n)), function(i) extract.clade(tre, node=i)$tip.label)
		for (i in which(sapply(spl, length) > (n/2))) {
			spl[[i]] <- setdiff(tre$tip.label, spl[[i]])
		}
		return(lapply(spl, sort))
	}
	map_splits <- function(query, reference) {
		m <- sapply(query, function(q, r) sapply(r, identical, y=q), r=reference)
		return(apply(m, 2, function(x) ifelse(any(x), which(x), NA)))
	}
	fun <- match.fun(method)
	class(trees) <- "multiPhylo"
	target <- splits(phy)
	sample <- lapply(trees, splits)
	matches <- sapply(seq_along(sample), function(i) map_splits(target, reference=sample[[i]]))
	edgelen <- sapply(seq_along(trees), function(i) trees[[i]]$edge.length[matches[,i]])
	phy$edge.length <- apply(edgelen, 1, fun, na.rm=TRUE)
	return(phy)
}


### FUNCTION
# prob_clades
### DESCRIPTION
# estimates posterior probability of clades in a tree from their frequency in posterior sample of trees
### ARGUMENTS
# phy: tree (class 'phylo') whose node heights are to be estimated
# psample: posterior sample of trees (class 'multiPhylo' or list of 'phylo' objects)
### VALUE
# numeric vector of posterior probabilities for clades ordered according to node numbers of their ancestors
prob_clades <- function(phy, psample) {
	clades <- ape::prop.part(phy)
	clades <- lapply(clades, function(x) attr(clades, "labels")[x])
	clades <- lapply(clades, sort)
	part <- ape::prop.part(psample)
	pppart <- attr(part, "number") / length(psample)
	part <- lapply(part, function(x) attr(part, "labels")[x])
	part <- lapply(part, sort)
	ord <- sapply(seq_along(clades),  function(i) which(sapply(part, identical, y=clades[[i]])))
	pp <- pppart[ord]
	return(pp)
}


### FUNCTION
# HPD
### DESCRIPTION
# estimates highest posterior density interval
### ARGUMENTS
# x: posterior sample of a model paramater
# p: the proportion of posterior density be included in the interval (defining 95% HPD by default, but possibly also 50% HPD and so on)
### VALUE
# numeric vector of posterior probabilities for clades ordered according to node numbers of their ancestors
### DETAILS
# "The 95% HPD stands for highest posterior density interval and represents the most compact interval on the selected parameter that contains 95% of the posterior probability. It can be loosely thought of as a Bayesian analog to a confidence interval." (Drummond & Bouckaert 2015, p. 89)

HPD <- function(x, p=0.95) {
	x <- na.rm(x)
	n <- round(length(x) * p)
	w <- seq(1, length(x) - n)
	x1 <- unname(sort(x))
	x2 <- rev(x1)
	int1 <- unname(rbind(x1[w + n], x1[w]))
	hpd1 <- sort(int1[, which.min(diff(int1))])
	int2 <- unname(rbind(x2[w + n], x2[w]))
	hpd2 <- sort(int2[, which.min(diff(int2))])
	hpd <- list(hpd1, hpd2)
	hpd <- hpd[[which.min(sapply(hpd, diff))]]
	return(hpd)
}


### FUNCTION
# discard_burnin
### DESCRIPTION
# discards a proportion of trees corresponding to a burnin of MCMC
### ARGUMENTS
# trees: object of class 'multiPhylo' or list of 'phylo' objects
# burnin: burn-in proportion of the posterior sample specified as a porportion, i.e. interval [0,1), or a percentage, interval [0, 100)
### VALUE
# the object 'trees' with the burnin part discarded

discard_burnin <- function(trees, burnin) {
	if (burnin == 0) {
		return(trees)
	} else if (burnin < 1) {
		return(trees[-seq(ceiling(burnin * length(trees)))])
	} else if (burnin < 100) {
		return(trees[-seq(ceiling(0.01 * burnin * length(trees)))])
	} else if (burnin < 0 | burnin > 100) {
		return(NA)
	}
}



### FUNCTION
# combine_runs
### DESCRIPTION
# combines posterior samples from independent runs of the analysis
### ARGUMENTS
# trees: list of 'multiPhylo' objects or lists of 'phylo' objects
# burnin: burn-in proportion of the posterior sample specified as a porportion, i.e. interval [0,1), or a percentage, interval [0, 100)
### VALUE
# the object of class 'multiPhylo' with posterior samples combined (possibly after discarding their burnin parts)

combine_runs <- function(trees, burnin=0) {
	if (burnin > 0) {
		trees <- lapply(trees, discard_burnin, burnin=burnin)
	}
	trees <- Reduce(c, trees)
	class(trees) <- "multiPhylo"
	return(trees)
}




### FUNCTION
# root_and_drop
### DESCRIPTION
# roots a tree or a series of them using specified outgroups and then removes the outgroups from the tree(s)
### ARGUMENTS
# phy: 'multiPhylo' object or a list of 'phylo' objects
# outgroup: character vector listing outgroup sequences or their unique identifier(s), e.g., a character string "outgroup"
### VALUE
# the object of the same class as the input 'phy'

root_and_drop <- function(phy, outgroup) {
	rootanddrop <- function (phy, outgroup) {
		return(ape::drop.tip(ape::root(phy, outgroup=outgroup, resolve.root=TRUE), tip=outgroup))
	}
	multi <- inherits(phy, "multiPhylo")
	if (multi | inherits(phy[[1]], "phylo")) {
		outgroup <- phy[[1]]$tip.label[sapply(outgroup, grep, x=phy[[1]]$tip.label)]
		phy <- lapply(phy, rootanddrop, outgroup=outgroup)
	} else {
		outgroup <- phy$tip.label[sapply(outgroup, grep, x=phy$tip.label)]
		phy <- rootanddrop(phy, outgroup=outgroup)
	}
	if (multi) {
		class(phy) <- "multiPhylo"
	}
	return(phy)
}


### FUNCTION
# remove_fossils
### DESCRIPTION
# removes fossil species identified as having terminal branches not reaching the maximum root-to-tip distance
### ARGUMENTS
# phy: 'multiPhylo' object or a list of 'phylo' objects
# tol: the maximum allowed difference from the maximum root-to-tip distance
### VALUE
# the object of the same class as the input 'phy'

remove_fossils <- function(phy, tol=1e-5) {
	rmfossils <- function (phy, tol) {
		n <- ape::Ntip(phy)
		dst <- ape::dist.nodes(phy)[n+1, 1:n]
		fossils <- which((max(dst) - dst) > tol)
		return(ape::drop.tip(phy, tip=fossils))
	}
	multi <- inherits(phy, "multiPhylo")
	if (multi | inherits(phy[[1]], "phylo")) {
		phy <- lapply(phy, rmfossils, tol=tol)
	} else {
		phy <- rmfossils(phy, tol=tol)
	}
	if (multi) {
		class(phy) <- "multiPhylo"
	}
	return(phy)
}



### REFERENCES
# Bouckaert R, Heled J, Kühnert D, Vaughan T, Wu CH, Xie D, Suchard MA, Rambaut A, Drummond AJ (2014) BEAST 2: a software platform for Bayesian evolutionary analysis. PLoS Computational Biology, 10: e1003537. https://doi.org/10.1371/journal.pcbi.1003537 
# Drummond AJ, Bouckaert R (2015) Bayesian evolutionary analysis with BEAST, Bayesian Evolutionary Analysis with BEAST. Cambridge University Press. https://doi.org/10.1017/CBO9781139095112
# Schliep KP (2011) phangorn: phylogenetic analysis in R. Bioinformatics, 27: 592-593. https://doi.org/10.1093/bioinformatics/btq706
# Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526-528. https://doi.org/10.1093/bioinformatics/bty633
