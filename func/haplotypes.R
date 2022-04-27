### haplotypes
# seqlen
### DESCRIPTION
# extracts unique haplotypes from an alignment
### ARGUMENTS
# input: name of the file (fasta or phylip interleaved format) or sequence alignment in 'DNAbin' or 'matrix' format
# fasta: optional, the name of fasta file with an alignement of haplotypes
# txt: optional, the name of txt file with information about assignment of sequences to haplotypes
# return: whether to return the output as an R object (a list)
# outgroup: character vector listing outgroup sequences or their unique identifier(s), e.g., a character string "outgroup"
# overlap: a proportion of overlap between sequence and haplotype, required to consider whether or not the sequence is assigned to the haplotype
### VALUE
# if return == TRUE, a list with components 'haplotypes' (the alignement of haplotypes), 'assign' (the assignment of sequences to haplotypes) & 'no' (the number of unique haplotypes found)

haplotypes <- function(input, fasta=NULL, txt=NULL, return=TRUE, outgroup=NULL, overlap=1) {
	if (inherits(input, "DNAbin")) {
		seq <- input <- toupper(as.character(input))
	} else if (inherits(input, "matrix")) {
		seq <- input <- toupper(input)
	} else {
		if (substr(readLines(input, n=1), 1, 1) == ">") {
			seq <- input <- toupper(ape::read.dna(input, format="fasta", as.character=TRUE, as.matrix=TRUE))
		} else {
			seq <- input <- toupper(ape::read.dna(input, as.character=TRUE, as.matrix=TRUE))
		}	
	}
	if (!is.null(outgroup)) {
		o <- sort(sapply(outgroup, grep, x=rownames(seq)))
		outgr <- seq[o,,drop=FALSE]
		seq <- seq[-o,,drop=FALSE]
		if (length(outgr) == 0) {
			outgroup <- NULL
			warning("alignment doesn't contain any of the specified outgroups")
		} else {
			outgroup <- outgr
		}
	}
	seq <- matrix(ifelse(seq %in% c("A","C","T","G"), seq, NA), nrow(seq), ncol(seq), dimnames=dimnames(seq))
	seq <- seq[order(rowSums(!is.na(seq)), decreasing=TRUE),]
	cvr <- rep(1, nrow(seq)) %*% t(rowSums(!is.na(seq)))
	ovr <- (!is.na(seq)) %*% t(!is.na(seq))
	ovr <- ovr / cvr
	dst <- ape::dist.dna(x=ape::as.DNAbin(seq), model="raw", pairwise.deletion=TRUE, as.matrix=TRUE)
	if (all(dst[upper.tri(dst)] > 0)) {
		haps <- setNames(vector("list", nrow(dst)), rownames(dst))
	} else {
		haps <- as.list(apply(dst == 0 & ovr >= overlap & upper.tri(dst), 1, which))
		haps <- lapply(haps, names)
		for (i in seq_along(haps)[-1]) {
			haps[[i]] <- setdiff(haps[[i]], unlist(haps[1:(i-1)]))
		}
		haps <- haps[!names(haps) %in% unlist(haps)]
	}
	haps <- Map(c, names(haps), haps)
	seq <- input[rownames(input) %in% sapply(haps,"[",1),,drop=FALSE]
	seq[is.na(seq)] <- "-"
	if (!is.null(outgroup)) {
		seq <- rbind(seq, outgroup)
	}
	if (!is.null(fasta) & is.null(txt)) {
		txt <- sub("\\.[[:alpha:]]+$", ".txt", fasta)
	}
	if (!is.null(txt) & is.null(fasta)) {
		fasta <- sub("\\.[[:alpha:]]+$", ".fasta", txt)
	}
	if (!is.null(fasta)) {
		ape::write.dna(x=seq, file=fasta, format="fasta")
	}
	if (!is.null(txt)) {
		writeLines(sapply(haps, paste, collapse=" "), txt)
	}
	if (isTRUE(return)) {
		return(list(haplotypes=seq, assign=haps, no=nrow(seq)))
	}
}



### FUNCTION
# seqlen
### DESCRIPTION
# calculates sequence lengths in the input alignment
### ARGUMENTS
# input: name of the file (fasta or phylip interleaved format) or sequence alignment in 'DNAbin' or 'matrix' format
### VALUE
# vector of sequence lengths
seqlen <- function(input) {
	if (inherits(input, "DNAbin")) {
		seq <- toupper(as.character(input))
	} else if (inherits(input, "matrix")) {
		seq <- toupper(input)
	} else {
		if (substr(readLines(input, n=1), 1, 1) == ">") {
			seq <- toupper(ape::read.dna(input, format="fasta", as.character=TRUE, as.matrix=TRUE))
		} else {
			seq <- toupper(ape::read.dna(input, as.character=TRUE, as.matrix=TRUE))
		}	
	}
	seq <- ifelse(seq == "-" | seq == "N" | seq == "n" | seq == "?", "-", seq)
	return(rowSums(seq != "-"))
}



### FUNCTION
# divergentseq
### DESCRIPTION
# selects given number of the most divergent sequences
### ARGUMENTS
# input: name of the file (fasta or phylip interleaved format) or sequence alignment in 'DNAbin' or 'matrix' format
# n: number of sequences required
# outgroup: character vector listing outgroup sequences or their unique identifier(s), e.g., a character string "outgroup"
# model: nucleotide substitution model used to calculate pairwise distences between sequences (deafult is 'raw', i.e. no correction for multiple mutations)
### VALUE
# an alignment (in 'DNAbin' or 'matrix' format) with the specificied no. of the most divergent sequences

divergentseq <- function(input, n, outgroup=NULL, model="raw") {
	if (inherits(input, "DNAbin")) {
		seq <- input <- toupper(as.character(input))
	} else if (inherits(input, "matrix")) {
		seq <- input <- toupper(input)
	} else {
		if (substr(readLines(input, n=1), 1, 1) == ">") {
			seq <- input <- toupper(ape::read.dna(input, format="fasta", as.character=TRUE, as.matrix=TRUE))
		} else {
			seq <- input <- toupper(ape::read.dna(input, as.character=TRUE, as.matrix=TRUE))
		}	
	}
	if (!is.null(outgroup)) {
		o <- sort(sapply(outgroup, grep, x=rownames(seq)))
		outgr <- seq[o,,drop=FALSE]
		seq <- seq[-o,,drop=FALSE]
		if (length(outgr) == 0) {
			outgroup <- NULL
			warning("alignment doesn't contain any of the specified outgroups")
		} else {
			outgroup <- outgr
		}
	}
	seq <- seq[order(seqlen(seq)),]
	dst <- ape::dist.dna(ape::as.DNAbin(seq), pairwise.deletion=TRUE, model=model, as.matrix=TRUE)
	ord <- which(lower.tri(dst), arr.ind=TRUE)[order(dst[lower.tri(dst)]),]
	excl <- minord <- min(ord[1,])
	while (length(excl) < (nrow(seq) - n)) {
		ord <- ord[!apply(ord == minord, 1, any),]
		minord <- min(ord[1,])
		excl <- c(excl, minord)
	}
	seq <- toupper(seq[-excl,])
	seq <- input[rownames(input) %in% c(rownames(seq), rownames(outgroup)),]
	return(seq)
}


