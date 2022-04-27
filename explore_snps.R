### TOOLS
library(ape)
library(dplyr)
source("func/snps_utils.R")
source("func/concatenate.R")
source("func/oneperind.R")
source("func/elbow.R")
source("func/writeSTRUCTURE.R")
source("func/branchcutting.R")
read.fasta <- function(file) toupper(ape::read.dna(file, format="fasta", as.matrix=TRUE, as.character=TRUE))
write.fasta <- function(seq, file) ape::write.dna(toupper(as.matrix(seq)), file, format="fasta")
write.delim <- function(x, file, row.names=FALSE) write.table(x, file, row.names=row.names, quote=FALSE, sep="\t")



### SNPS
# loading of pre-selected loci, calling of biallelic SNPs (one per locus)
# preparing input files for different analyses
# loading of information about classification of individuals (missing in purely exploratory analyses)

load("data/geo_loci_subset.R")

geo_usnps <- get_snps(geo_loci_subset, nalleles=2, type="pis", unique=TRUE, as.matrix=FALSE)
geo_binary <- make_binary_snps(concatenate(geo_usnps), center=TRUE, type="snps", onerowperind=TRUE, allele="_A[[:digit:]]$")
geo_concat <- concatenate(oneperind(geo_usnps, allele="_A[[:digit:]]$", remove=TRUE))
write.fasta(geo_concat, "data/georychus_unsps_concat.fasta")
geo_strdat <- write_str_data(loci=geo_usnps, file="data/georychus_unsps.str", allele="_A[[:digit:]]$", return=TRUE)

geo_info <- read.delim("data/georychus_info.txt", stringsAsFactors=FALSE)
geo_palette <- c(Cape="#e6194b", Struisbaai="#f58231", Oudshoorn="#3cb44b", KwaZuluNatal="#4363d8", Mpumalanga="#46f0f0")


### PCA
# PCA on counts of biallelic SNPs at different loci as input variables
# it is performed without scaling as all variables are counts (and hence mutually comparable)
# and without centering as they are pre-centered so the heterozygote (natural centers) ~ 0
# elbow method of Salvador & Chan (2004) is used to determine number of PCs to be retained (k)
# PS scores are displayed for PC1-2, for all pairs of PC1-k and then used for UPGMA clustering

geo_counts <- count_binary_snps(geo_binary)
geo_pca <- prcomp(geo_counts, center=FALSE, scale.=FALSE)
geo_pcs <- geo_pca$x
geo_bg <- geo_palette[geo_info$Lineage[match(rownames(geo_pcs), geo_info$ddRAD)]]
k <- elbow(geo_pca$sdev^2)$best

par(mai=c(1.02,1.02,0.42,0.42))
plot(geo_pcs, pch=21, cex=2, bg=geo_bg, cex.lab=1.5, cex.axis=1.25)

pairs(geo_pcs[,1:k], pch=21, cex=2, bg=geo_bg, lwd=0.5)
clust <- hclust(dist(geo_pcs[,1:k], method="manhattan"), method="average")
plot(clust)


### CONCATENATED ML TREE
# assumes 'iqtree' binary in 'bin' folder
# ML tree is estimated and partitioned into OTUs using branchcutting method
# ML tree is displayed with tips colored according to OTUs

iqtree <- ifelse(platform == "unix", "./bin/iqtree2", "bin/iqtree.exe")
system2(iqtree, args=c("-s", "data/georychus_unsps_concat.fasta", "-m", "MFP", "-B", "1000", "-T", "AUTO", "--prefix", "results/georychus_unsps_concat"))

mltree <- ape::read.tree("results/georychus_unsps_concat.treefile")
bcut <- branchcutting(mltree)
otus <- write.bcut(bcut, "results/georychus_bcut_otus.txt", return=TRUE)

geo_bg <- geo_palette[geo_info$Lineage[match(mltree$tip.label, geo_info$ddRAD)]]
#geo_bg <- otus$OTU[match(mltree$tip.label, otus$ID)]
plot(mltree, type="unrooted", lab4ut="axial", no.margin=TRUE, label.offset=0.01, cex=0.7)
ape::tiplabels(pch=21, col=1, bg=geo_bg, cex=1.5)



### STRUCTURE
# parses STRUCTURE control files
# estimates admixture coefficients for the specified K (here 5)
# loads the results and displays them in the form of barplot with color opaqueness indicating confidence about the estimates (argument 'INTERVAL=TRUE')

K <- 5
geo_strpar <- 
data.frame(
Arg=c("infile", "outfile", "MAXPOPS", "NUMINDS", "NUMLOCI", "ONEROWPERIND", "LABEL", "POPDATA", "POPFLAG", "LOCDATA", "MAPDISTANCES", "NOADMIX", "LINKAGE", "USEPOPINFO", "LOCPRIOR", "LOCISPOP", "ANCESTDIST", "ANCESTPINT","BURNIN","NUMREPS"),
Value=c("data/georychus_unsps.str", "results/georychus_unsps.res", K, nrow(geo_strdat), ncol(geo_strdat)/2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.90, "30000", "100000"))

parse_str_params(df=geo_strpar, templates=c("data/mainparams", "data/extraparams"), outputs=c("data/georychus_mainpar.txt", "data/georychus_extrapar.txt"))

structure <- ifelse(platform == "unix", "./bin/structure", "bin/structure.exe")
system2(structure, args=c("-m", "data/georychus_mainpar.txt", "-e", "data/georychus_extrapar.txt"))

geo_str_res <- read_str_output("results/georychus_unsps.res_f")
cluster <- apply(geo_str_res[,seq(attr(geo_str_res,"K"))], 1, which.max)
mscspec <- geo_info$Lineage[match(names(cluster), geo_info$ddRAD)]
matching <- sort(apply(table(mscspec, cluster), 1, which.max))
plot_str_output(geo_str_res, interval=TRUE, palette=geo_palette[names(matching)])


