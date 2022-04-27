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

