read_raxml_partitions <- function(file) {
	p <- Filter(function(x) nchar(x) > 0, readLines(file))
	p <- gsub("=", "-", gsub("^DNA[[:alpha:]]*,|\\s", "", p))
	p <- do.call(rbind, strsplit(p, "-"))
	return(data.frame(start=as.numeric(p[,2]), end=as.numeric(p[,3]), row.names=p[,1]))
}

write_raxml_partitions <- function(p, file) {
	if (ncol(p) == 2) {
		p$part <- "123"
	}
	p[,3] <- ifelse(is.na(p[,3]), "123", p[,3])
	p[,3] <- ifelse(nchar(p[,3]) == 3, "123", p[,3])
	codon <- strsplit(p[,3], "[[:punct:]]")
	times <- c(1,3)[(sapply(codon, length) > 1) + 1]
	p <- split(p[rep(seq(nrow(p)), times),1:2], rep(rownames(p), times))[rownames(p)]
	for (i in seq_along(p)) {
		if (length(codon[[i]]) == 1) {
			p[[i]] <- paste("DNA,", paste(names(p)[i], "=", paste0(p[[i]][,1], "-", p[[i]][,2])))
		} else {
			partnam <- paste0(names(p)[i], codon[[i]])
			codon[[i]] <- strsplit(codon[[i]], "")
			p[[i]][,1] <- p[[i]][,1] + (as.numeric(unlist(codon[[i]])) - 1)
			p[[i]] <- paste0(p[[i]][,1], "-", p[[i]][,2], "\\3")	
			p[[i]] <- split(p[[i]], rep(seq_along(codon[[i]]), sapply(codon[[i]], length)))
			p[[i]] <- sapply(p[[i]], paste, collapse=", ")
			p[[i]] <- paste("DNA,", partnam, "=", p[[i]])
		}
	}
	p <- unname(unlist(p))
	writeLines(p, file)
}
