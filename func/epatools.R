### FUNCTION
# read_jplace
### DESCRIPTION
# reads in an EPA output in .jplace format
### ARGUMENTS
# jp: .jplace output of EPA/pplacer
### VALUE
# object of class 'jplace', which is a list with components:
	# tree: is an unrooted tree EPA was operating on (in 'phylo' format)
	# placements: is a data.frame with EPA output, i.e., list of possible placements of query sequences into the tree and their supports (likelihoods or posterior probabilities)
	# software: indicates software used (currently, it assumes either 'pplacer' or 'raxml')

read_jplace <- function(jp) {
# function
	getstring <- function(pattern, text) regmatches(text, regexpr(pattern, text)) 
# input
	jp <- gsub("\\s", "", readLines(jp))
# tree
	tree <- paste0("(", gsub("^.*\"\\(|;.*$", "", jp[grep("\"tree\"", jp)]), ";")
	edge <- gregexpr("\\{[[:digit:]]+}", tree)
	edgelabel <- Map(paste0, regmatches(tree, edge), ":")
	regmatches(tree, edge) <- ""
	node <- gregexpr(":", tree)
	regmatches(tree, node) <- edgelabel
	tree <- ape::read.tree(text=tree)
	tiplabel <- regmatches(tree$tip.label, gregexpr("\\{\\d+}", tree$tip.label))
	tiplabel <- gsub("[\\{}]", "", unlist(tiplabel))
	nodelabel <- regmatches(tree$node.label, gregexpr("\\{\\d+}", tree$node.label))
	nodelabel <- lapply(nodelabel, function(x) gsub("[\\{}]", "", x))
	nodelabel <- unlist(ifelse(sapply(nodelabel, length) == 0, NA, nodelabel))
	edgelabel <- c(tiplabel, nodelabel)
	tree$tip.label <- gsub("\\{[[:digit:]]+}", "", tree$tip.label)
	tree$node.label <- NULL
	tree$edge.label <- edgelabel[tree$edge[,2]]
# placements
	place <- jp[grep("\"p\":", jp)]
	place <- setNames(gsub("^\\[|]$", "", regmatches(place, regexpr("\\[\\[.*]]", place))),
		gsub("\\[\"|\"]", "", regmatches(place, regexpr("\\[\".*\"]", place))))
	place <- lapply(place, function(p) unlist(regmatches(p, gregexpr("\\[.*]{1}?", p))))
	place <- lapply(place, function(p) gsub("^\\[|]$", "", unlist(p)))
	place <- lapply(place, function(p) do.call(rbind, lapply(strsplit(p, ","), as.numeric)))
	fields <- jp[grep("fields", jp):length(jp)]
	fields <- paste(fields[Reduce(":", grep("\\[|]", fields))], collapse="")
	fields <- unlist(strsplit(gsub("^.*\\[|].*$|\"", "", fields), ","))
	place <- setNames(data.frame(rep(names(place), sapply(place, nrow)), do.call(rbind, place)), c("id", fields))
# software
	soft <- jp[grep("\"invocation\"", jp)]
	soft <- unlist(lapply(c("raxml","pplacer"), getstring, text=tolower(soft)))
# output
	result <- list(tree=tree, placements=place, software=soft)
	attr(result, "classified") <- FALSE
	class(result) <- "jplace"
	return(result)
}


### FUNCTION
# classify_jplace
### DESCRIPTION
# classifies edges with phylogenetic placements into clades / bipartitions
### ARGUMENTS
# jp: object of class 'jplace' (= output of read_jplace), i.e. list with components:
	# tree: object of class 'phylo' amended by edge.label component
	# placements: data.frame
	# software: name of the software used (currently pplacer or raxml)
# clade: data frame with components 'tip' and 'clade' defining classification of tips to (possibly nested) clades / bipartitions (if the names differ, the first column is considered as 'tip' and the second as 'clade')
# ancestral: whether to classify ancestral branches (only if the tree is rooted and defined clades are comprehensive and non-overlapping, otherwise switched to FALSE and issues a warning)
### VALUE
# an amended 'jplace' object with additional components 'clade' and 'position' in the 'placements' data frame indicating classification of the branches
# 


classify_jplace <- function(jp, clade, ancestral=TRUE) {
	getpaths <- function(phy, from, to) lapply(to, function(x) ape::nodepath(phy, from=from, to=x))
	find_common_branch <- function(tip, phy) ifelse(length(tip) == 1, match(tip, phy$tip.label), ape::getMRCA(phy, tip))
	tree <- jp$tree
	place <- jp$placements
	names(clade) <- tolower(names(clade))
	if (any(!c("tip", "clade") %in% names(clade))) {
		names(clade)[1:2] <- c("tip", "clade")	
	}
	clade <- clade[clade$tip %in% jp$tree$tip.label,]
	clade <- split(clade$tip, clade$clade)
	anc <- lapply(clade, find_common_branch, phy=tree)
	off <- lapply(clade, match, table=tree$tip.label)
	paths <- Map(getpaths, list(tree), anc, off)
	paths <- lapply(lapply(lapply(paths, unlist), unique), sort)
	paths <- Map(setdiff, paths, anc)
	crowns <- setNames(lapply(paths, function(x) tree$edge.label[match(x, tree$edge[,2])]), names(clade))
	crown_lengths <- sapply(crowns, length)
	crowns <- data.frame(edge=unlist(crowns), clade=rep(names(crowns), crown_lengths), position=rep("crown", sum(crown_lengths)), row.names=NULL)	
	stems <- setNames(tree$edge.label[match(unlist(anc), tree$edge[,2])], names(clade))
	stems <- data.frame(edge=stems, clade=names(stems), position="stem", row.names=NULL)
	outgroup_tips <- as.numeric(sapply(tree$outgroup, grep, x=tree$tip.label))
	overlapping <- any(duplicated(unlist(paths)))
	partial <- any(!setdiff(seq(ape::Ntip(tree)), outgroup_tips) %in% unlist(off))
	unrooted <- !ape::is.rooted(tree)
	if (ancestral & (overlapping | unrooted | partial)) {
		ancestral <- FALSE
		message <- c("defined clades are overlapping", "the classification is not comprehensive", "the tree is not rooted")
		message <- message[c(overlapping, partial, unrooted)]
		if (length(message) > 1) {
			message <- paste(paste(message[1:(length(message) - 1)], collapse=", "), message[length(message)], sep=" and ")	
		}
		message <- paste0(toupper(substr(message, 1, 1)), substr(message, 2, nchar(message)), ",")
		warning(paste(message, "'ancestral' argument was set to FALSE."))
	}
	if (ancestral) {
		outgroup_paths <- getpaths(tree, from=ape::Ntip(tree)+1, to=outgroup_tips)
		root_neighbors <- tree$edge[tree$edge[,1] == ape::Ntip(tree) + 1,2]
		excluded <- c(unlist(anc), unlist(paths), unlist(outgroup_paths), root_neighbors)
		mrca <- setdiff(seq(ape::Nnode(tree))[-1] + ape::Ntip(tree), excluded)
		ancestral_tips <- lapply(mrca, function(x) ape::extract.clade(tree, node=x)$tip.label)
		ancestral_tips <- lapply(ancestral_tips, match, table=tree$tip.label)
		ancestral_paths <- Map(getpaths, list(tree), as.list(mrca), ancestral_tips)
		ancestral_paths <- lapply(lapply(lapply(ancestral_paths, unlist), unique), sort)
		superspecies <- lapply(ancestral_paths, function(x) sapply(anc, "%in%", x))
		superspecies <- lapply(superspecies, function(x) names(x[x]))
		superspecies <- sapply(superspecies, paste, collapse="_")
		ancestors <- setNames(tree$edge.label[match(mrca, tree$edge[,2])], superspecies)
		ancestors <- data.frame(edge=ancestors, clade=names(ancestors), position="ancestral", row.names=NULL)
	} else {
		ancestors <- data.frame(edge=character(0), clade=character(0), position=character(0))
	}
	edges <- do.call(rbind, list(crowns, stems, ancestors))
	ii <- match(place$edge_num, edges$edge)
	place$clade <- edges$clade[ii]
	place$position <- edges$position[ii]
	jp$placements <- place
	attr(jp, "classified") <- TRUE
	class(jp) <- "jplace"
	return(jp)
}


### FUNCTION
# root_jplace
### DESCRIPTION
# roots the tree & adjusts edge labelling accordingly
# ARGUMENTS
# jp: object of class "jplace" (= output of read_jplace)
# outgroup: character vector listing outgroup sequences or their unique identifier(s), e.g., a character string "outgroup"
### VALUE
# a modified object of class 'jplace', i.e., its 'tree' is rooted and 'placements$edge_num' modified

root_jplace <- function(jp, outgroup) {
	if (attr(jp, "classified") == TRUE) {
		stop("classified 'jplace' objects are not assumed to be rooted")
	}
	tree <- jp$tree
	match_tips <- match(seq(ape::Ntip(tree)), tree$edge[,2])
	match_nodes <- match(seq(ape::Nnode(tree)) + ape::Ntip(tree), tree$edge[,2])
	tree$tip.label <- paste0(tree$tip.label, "{", tree$edge.label[match_tips], "}")
	tree$node.label <- paste0("{", tree$edge.label[match_nodes], "}")
	outgroup_tips <- tree$tip.label[sapply(outgroup, grep, x=tree$tip.label)]
	tree <- ape::root(tree, outgroup=outgroup_tips, resolve.root=TRUE)	
	tree$edge.label <- character(ape::Nedge(tree))
	match_tips <- match(seq(ape::Ntip(tree)), tree$edge[,2])
	match_nodes <- match(seq(ape::Nnode(tree)) + ape::Ntip(tree), tree$edge[,2])
	tree$edge.label[match_tips] <- tree$tip.label
	tree$edge.label[na.omit(match_nodes)] <- tree$node.label[!is.na(match_nodes)]
	tree$edge.label <- gsub("^.*\\{|}", "", tree$edge.label)
	tree$tip.label <- gsub("\\{[[:digit:]]+}$", "", tree$tip.label)
	tree$node.label <- NULL
	root_edge <- which(tree$edge[,1] == ape::Ntip(tree) + 1)
	root_edge_label <- na.omit(suppressWarnings(as.numeric(tree$edge.label[root_edge])))
	tree$edge.label[root_edge] <- NA
	tree$outgroup <- outgroup
	jp$tree <- tree
	jp$placements$edge_num[jp$placements$edge_num %in% root_edge_label] <- NA
	attr(jp, "classified") <- FALSE
	class(jp) <- "jplace"
	return(jp)
}


### FUNCTION
# classify_sequences
### DESCRIPTION
# classifies sequences to clades / bipartitions
# ARGUMENTS
# jp: object of class 'jplace' pre-treated by 'classify_jplace' (i.e. containing classification to clades / bipartitions)
# keep: no. of placements per clade to keep (keep=1 - ML placement, keep=NULL means all), corresponds to 'epa-­keep-­placements' in RAxML
# probthr: likelihood weight threshold, corresponds to 'epa­-prob­-threshold' in RAxML
# cumprobthr: cummulative likelihood weight threshold, corresponds to 'epa­-accumulated­-threshold' in RAxML (should be < 1 if there is not 'clade' variable in jp$placements)
# VALUE
# data frame with components:
	# id: label of query sequence
	# clade: label of clade / bipartition
	# position: crown or stem 
	# probability: probability of classification, i.e., sum of probabilities over the branches belonging to the clade / bipartition. In thos context, 'probability' means either posterior probability or relative likelihood weight, i.e., a quantity which sums up to unity over all branches of the tree.
	# n_branch: number of branches which contributed to the sum

classify_sequences <- function(jp, keep=NULL, probthr=0, cumprobthr=1) {
	if (attr(jp, "classified") == FALSE) {
		stop("the 'jplace' object was not pre-treated by 'classify_jplace'")
	}
	if (is.null(keep)) {
		keep <- ape::Nedge(jp$tree)
	}
	probability <- c(pplacer="post_prob", raxml="like_weight_ratio")[jp$software]
	place <- jp$placements
	place$label <- paste(place$clade, place$position, sep="_")
	place <- split(place, place$id)
	for (i in seq_along(place)) {
		iprob <- place[[i]][,probability]
		ikeep <- which(iprob >= probthr & cumsum(iprob) <= cumprobthr)
		if (length(ikeep) > 0) {
			ikeep <- seq(min(c(max(ikeep), keep)))
		} else {
			ikeep <- FALSE
		}
		place[[i]] <- place[[i]][ikeep,,drop=FALSE]
		if (nrow(place[[i]]) > 0) {
			place[[i]] <- split(place[[i]], place[[i]]$label)
			for (j in seq_along(place[[i]])) {
				prob <- data.frame(probability=sum(place[[i]][[j]][,probability]), n_branch=nrow(place[[i]][[j]]))
				place[[i]][[j]] <- cbind(place[[i]][[j]][1, c("id","clade","position"), drop=FALSE], prob)
			}
			place[[i]] <- do.call(rbind, place[[i]])
			place[[i]] <- place[[i]][order(place[[i]]$probability, decreasing=TRUE),,drop=FALSE]
		}
	}
	place <- do.call(rbind, place)
	rownames(place) <- NULL
	return(place)
}


### FUNCTION
# plot.jplace
### DESCRIPTION
# simple plot of the estimated placements, for query sequences (or possibly a subset of them) it highlights branches where they are probably placed
# ARGUMENTS
# jp: object of class "jplace" (= output of read_jplace)
# subset: vector of sequence names
# keep: how many first placements to keep for display
# prob: cummulative probability threshold up to which the placements are displayed
# type: type of tree, default is "radial" (branch lengths discraded!), but "unrooted" is also reasonable
# col: color to highlight placements
# bg: color of background branches
# file: file name, if specified the plot is save as a .pdf file of the given name
# ... other parameters of plot.phylo

plot.jplace <- function(jp, subset=NULL, keep, prob=1.00, type="radial", col=3, bg=1, file, ...) {
	probability <- c(pplacer="post_prob", raxml="like_weight_ratio")[jp$software]
	place <- jp$placements
	if (!probability %in% names(place)) {
		stop("edges of the tree in 'jplace' object are yet to be classified into clades")
	}
	if (!is.null(subset)) {
		place <- place[place$id %in% subset,,drop=FALSE] 		
	}
	place <- split(place, place$id)
	if (prob < 1) {
		place <- lapply(place, function(x) x[1:min(which(cumsum(x[,probability]) >= prob)),,drop=FALSE])
	}
	if (!missing(keep)) {
		place <- lapply(place, function(x) x[1:min(keep,nrow(x)),,drop=FALSE])		
	}
	place <- do.call(rbind, place)
	edgecol <- c(bg, col)[jp$tree$edge.label %in% place$edge_num + 1]
	if (!missing(file)) {
		pdf(file)
	}	
	plot(jp$tree, type=type, edge.color=edgecol, ...)
	if (!missing(file)) {
		dev.off()
	}	
}


### FUNCTION
# make_plotlabels
### DESCRIPTION
# makes labels indicating which query sequences are most probably placed on given branches
# ARGUMENTS
# jp: 'jplace' object with classified edges
# classification: output of 'classify_sequences'
# sep: a string separating labels of query sequences placed on the same branch
# VALUE
# character vector with concatenated id's of most probably placed query sequences; the vector elements are named by branch labels from 'jp$tree$edge.label'

make_plotlabels <- function(jp, classification, sep="\n") {
	best <- classification[match(unique(classification$id), classification$id),]
	probability <- c(pplacer="post_prob", raxml="like_weight_ratio")[jp$software]
	place <- jp$placements
	place <- split(place, place$id)
	best <- best[match(names(place),best$id),]
	for (i in seq_along(place)) {
		place[[i]] <- place[[i]][place[[i]]$clade == best$clade[i] & place[[i]]$position == best$position[i],,drop=FALSE]
		place[[i]] <- place[[i]][which.max(place[[i]][,probability]),,drop=FALSE]
	}
	place <- do.call(rbind, place)
	edges <- split(place$id, place$edge_num)
	edgelab <- sapply(edges, paste, collapse=sep)
	edgelab <- setNames(edgelab, match(names(edgelab), jp$tree$edge.label))
	return(edgelab)
}


### HELPER FUNCTIONS writing .fasta file and tab-delimited text file
write.fasta <- function(x, file) ape::write.dna(x, file, format="fasta")
write.delim <- function(x, file) write.table(x, file, quote=FALSE, row.names=FALSE, sep="\t")


### FUNCTION
# ml_classification
### DESCRIPTION
# finding maximum likelihood classification
### ARGUMENTS
# place: output of classify_sequences
### VALUE
# data frame with ML classification for each query sequence

ml_classification <- function(place) {
	mlplace <- place[match(unique(place$id), place$id),]
	mlplace$totalprob <- mlplace$probability
	for (i in seq(nrow(mlplace))) {
		matchid <- place$id == mlplace$id[i]
		matchclade <- place$clade == mlplace$clade[i]
		mlplace$totalprob[i] <- sum(place$probability[which(matchid & matchclade)])
	}
	return(mlplace)
}
