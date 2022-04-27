### NOTES & COMMENTS

# The pipeline:
# - selects a number of the most divergent unique haplotypes (longer than specified length)
# - infers maximum likelihood (ML) tree from the selected sequences (substitution model selection included)
# - delimits operational taxonomic units (OTUs)
# - estimates ML phylogenetic placement of the remaining (=query) sequences
# - classifies query sequences into the OTUs
# - displays OTU-colored ML tree
# - optionally displays OTU-colored map of sampling locations

# In its basic form it requires specifying just the name of input alignement and outgroup indicator,
# for the mapping and info file with geographical coordinates has to be supplied.
# The remaining sections (introduced by '###' header) can be run separately or in a series.
# Every section starts with automated specification of necessary file names,
# and then it calls either R function or third party software (IQTREE, mptp, RAxML).
# Also inputs and outputs are managed along the way
# (e.g., partition file is created, tab-delimited file is exported etc.)

# The worked out example is based on data from Uhrová et al. (2022)



### TOOLS
library(ape)
library(maps)
source("func/haplotypes.R")
source("func/raxml_partitions.R")
source("func/epatools.R")
source("func/read_ptp.R")
source("func/match_otus.R")
source("func/rootanddrop.R")
source("func/branchcutting.R")
source("func/plot_clades.R")
source("func/plot_phylomap.R")
read.fasta <- function(file) toupper(ape::read.dna(file, format="fasta", as.matrix=TRUE, as.character=TRUE))
write.fasta <- function(seq, file) ape::write.dna(toupper(as.matrix(seq)), file, format="fasta")
write.delim <- function(x, file, row.names=FALSE) write.table(x, file, row.names=row.names, quote=FALSE, sep="\t")



### SEQUENCE FILE NAME & OTHER PARAMETERS
# seqfile: the name of input alignment
# outgroup: names of outgroup sequences or their unique identifier(s), e.g., "outgroup"
# minlength: minimum length for sequences to be included in the tree, zero indicates including everything
# no: number of the most divergent haplotypes to be retained. 'Inf' means "retain all"
# platform: either 'unix' or 'windows'
# coding: whether to consider partitioning by codon position
# method: OTU-picking method, either "bcut" (in R, all platforms) or "mptp" (calling 'mptp', unix only)
# infofile: tab-delimited file with geographical coordinates, necessary for 'plotPhylogeo'

seqfile <- "data/heliophobius_cytb.fasta"
outgroup <- "Georychus"
minlength <- 700
no <- 150
platform <- .Platform$OS.type
coding <- TRUE
method <- "bcut"
infofile <- "data/heliophobius_info.txt"


### HAPLOTYPES
# data
cytb <- read.fasta(seqfile)
outgroups <- rownames(cytb)[grepl(outgroup, rownames(cytb))]
# file names
hapmap <- sub("\\..+$", "_haps.txt", seqfile)
hapseq <- sub("txt$", sub("^.+\\.", "", seqfile), hapmap)

# unique haplotypes
haplist <- haplotypes(seqfile, outgroup=outgroups)
writeLines(sapply(haplist$assign, paste, collapse=" "), hapmap)

# most divergent sequences
haps <- haplist$haplotypes[seqlen(haplist$haplotypes) >= minlength,,drop=FALSE]
if (no < (nrow(haps) - length(outgroups))) {
	haps <- divergentseq(input=haps, n=no, outgroup=outgroups, model="raw")
}
write.fasta(haps, hapseq)



### MAXIMUM LIKELIHOOD TREE (IQ-TREE)
# file names & parsing partition file
path <- regmatches(seqfile, regexpr("^.*\\/", seqfile))
if (length(path) == 0) path <- getwd()
hapseq <- paste0(sub("\\..*$", "", seqfile), "_haps.", sub("^.+\\.", "", seqfile))
prefix <- paste0(sub("\\..+$", "", seqfile), "_iqt")
partfile <- sub("\\..+$", "_part.txt", seqfile)
nbp <- ncol(read.fasta(hapseq))
part <- ifelse(isTRUE(coding), "1-2-3", "123")
write_raxml_partitions(data.frame(start=1, end=nbp, part=part, row.names="gene"), partfile)

# iq-tree
iqtree <- ifelse(platform == "unix", "./bin/iqtree2", "bin/iqtree.exe")
system2(iqtree, args=c("-s", hapseq, "-p", partfile, "-m", "TESTMERGE", "-mset", "HKY", "-mfreq", "F", "-mrate", "E,G", "-B", "1000", "-T", "4", "--prefix", prefix))

# post-processing of outputs
outputs <- list.files(sub("/$", "", path), pattern=sub(path, "", prefix, fixed=TRUE), full.names=TRUE)
best_scheme <- outputs[grep("\\.best_scheme$", outputs)]
treefile <- outputs[grep("\\.treefile$", outputs)]
invisible(file.remove(setdiff(outputs, c(best_scheme, treefile))))
invisible(file.rename(treefile, sub("_iqt", "", treefile)))
invisible(file.rename(best_scheme, partfile))
writeLines(gsub("DNAF", "DNA", readLines(partfile)), partfile)



### OTU PICKING (mPTP / branchcutting)
# file names
mptp_prefix <- paste0(sub("\\..+$", "", seqfile), "_mptp")
bcut_prefix <- paste0(sub("\\..+$", "", seqfile), "_bcut")
mptpfile <- paste0(mptp_prefix, ".txt")
bcutfile <- paste0(bcut_prefix, ".txt")
otufile <- sub("\\..+$", "_otus.txt", seqfile)

path <- regmatches(seqfile, regexpr("^.*\\/", seqfile))
if (length(path) == 0) path <- getwd()
treefile <- list.files(sub("/$", "", path), pattern="\\.treefile$", full.names=TRUE)
rtreefile <- sub("\\.treefile$", "_rooted.tre", treefile)
tree <- ape::read.tree(treefile)
rtree <- rootanddrop(tree, outgroup=outgroups)
ape::write.tree(rtree, rtreefile)
outgroups <- tree$tip.label[grepl(outgroup, tree$tip.label)]

# mptp (if platform == "unix") & branchcutting
if (platform == "unix") {
	system2("./bin/mptp", args=c("--ml", "--multi", "--tree_file", rtreefile, "--output_file", mptp_prefix))
} else if (method == "mptp") {
	method <- "bcut"
	warning("on a unix-type platform 'method' was switched to 'bcut'")
}
bcut <- branchcutting(rtree)
write.bcut(bcut, bcutfile, return=FALSE)

# post-processing of outputs
if (platform == "unix") {
	otus_mptp <- read_ptp(mptpfile)
	otus_bcut <- read.delim(bcutfile)
	otus <- match_otus(otus_mptp, otus_bcut, col.names=c("mptp", "bcut"))
} else {
	otus <- read.delim(bcutfile, col.names=c("ID", "bcut"))
}
write.delim(otus, otufile)



### PHYLOGENETIC PLACEMENTS (EPA)
# file names and copying
path <- regmatches(seqfile, regexpr("^.*\\/", seqfile))
if (length(path) == 0) path <- getwd()
prefix <- gsub("^.*\\/|\\..+$", "", seqfile)
treefile <- list.files(sub("/$", "", path), pattern="\\.treefile$", full.names=TRUE)
partfile <- sub("\\..+$", "_part.txt", seqfile)
lseqfile <- sub("^/", "", sub(path, "", seqfile))
ltreefile <- sub("^/", "", sub(path, "", treefile))
lpartfile <- sub("^/", "", sub(path, "", partfile))
jplacefile <- paste0(sub("/$", "", path), "/", prefix, ".jplace")
invisible(file.copy(seqfile, lseqfile))
invisible(file.copy(treefile, ltreefile))
invisible(file.copy(partfile, lpartfile))

# EPA (in RAxML)
raxml <- ifelse(platform == "unix", "./bin/raxmlHPC-SSE3", "./bin/raxmlHPC.exe")
system2(raxml, args=c("-f", "v", "-m", "GTRGAMMA", "--HKY85", "--epa-keep-placements=100", "--epa-prob-threshold=0.01", "-s", lseqfile, "-t", ltreefile, "-q", lpartfile, "-n", prefix))

# post-processing of outputs
outputs <- c(paste0(c(lseqfile, lpartfile), ".reduced"), list.files(pattern=paste0("RAxML.*\\.", prefix)))
invisible(file.rename(outputs[grep("\\.jplace$", outputs)], jplacefile))
toberemoved <- c(lseqfile, ltreefile, lpartfile)
nottoberemoved <- c(seqfile, sub("^/", "", sub(getwd(), "", treefile)), partfile)
invisible(suppressWarnings(file.remove(c(setdiff(toberemoved, nottoberemoved), outputs[-grep("\\.jplace$", outputs)]))))



### CLASSIFICATION
# file names
jplacefile <- sub("\\..+$", ".jplace", seqfile)
otufile <- sub("\\..+$", "_otus.txt", seqfile)

# read in .jplace file, root the tree & classify its branches & query sequences to OTUs
jplace <- read_jplace(jplacefile)
outgroups <- jplace$tree$tip.label[grepl(outgroup, jplace$tree$tip.label)]
jplace <- root_jplace(jplace, outgroup=outgroups)
otus <- read.delim(otufile)
otus <- otus[,c(names(otus)[1], method)]
jplace_otus <- classify_jplace(jplace, otus, ancestral=TRUE)
class_otus <- classify_sequences(jplace_otus)
class_otus_ml <- ml_classification(class_otus)

# post-processing outputs
otus <- data.frame(otus, position="crown", probability=1, query=FALSE, stringsAsFactors=FALSE)
query <- data.frame(class_otus_ml[,c("id","clade","position","probability")], query=TRUE, stringsAsFactors=FALSE)
names(query)[1:4] <- c("ID", names(otus)[2], "position", "probability")
otus <- rbind(otus, query)
write.delim(otus, otufile)



### PLOT TREE
# file names
rtreefile <- sub("\\..+$", "_rooted.tre", seqfile)
otufile <- sub("\\..+$", "_otus.txt", seqfile)

# data & plot
mltree <- ape::read.tree(rtreefile)
otus <- read.delim(otufile)
plot_clades(mltree, otus)


### PLOT MAP
# file names
otufile <- sub("\\..+$", "_otus.txt", seqfile)

# data & map
otus <- read.delim(otufile)
info <- read.delim(infofile)
plot_phylomap(info, otus)

