#!/usr/bin/env Rscript

########################################
## hacked script to run dupRadar

## http://bioconductor.org/packages/release/bioc/vignettes/dupRadar/inst/doc/dupRadar.html
## https://bioconductor.org/packages/devel/bioc/manuals/Rsubread/man/Rsubread.pdf

## input:
## - _duplicate marked_ bam file
## - gtf file
## - parameters for duplication counting routine:
##   stranded, paired, outdir, threads.

## make sure the input bam has marked duplicates (using picard markduplicates or biobambam2)
#picard.sh I=bt2.output.bam.q10_filter.sorted.bam O=dupmarked.bam M=marked.txt


########################################
## set path to user's R libraries
libLocation <- "/pub46/willr/R/x86_64-unknown-linux-gnu-library/3.2"


########################################
## load dupRadar
.libPaths( libLocation )
library(dupRadar)


########################################
## parse CLI input

args   <- commandArgs(TRUE)

# input bam file (duplicates must be marked)
bam <- gsub("bam=","",args[1])

# GFF annotation file
gtf <- gsub("gtf=","",args[2])

# stranded library? no|yes|reverse
stranded <- gsub("stranded=","",args[3])

# paired end experiment?
paired   <- gsub("paired=","",args[4])

# output directory
outdir   <- gsub("outdir=","",args[5])

# number of threads
threads  <- as.integer(gsub("threads=","",args[6]))

# command checks:
if(!file.exists(bam)) {
  stop(paste("File",bam,"does NOT exist"))
}

if(!file.exists(gtf)) {
  stop(paste("File",gtf,"does NOT exist"))
}

if(!file.exists(outdir)) {
  stop(paste("Dir",outdir,"does NOT exist"))
}

if(is.na(stranded) | !(grepl("no|yes|reverse",stranded))) {
  stop("Stranded has to be no|yes|reverse")
}

if(is.na(paired) | !(grepl("no|yes",paired))) {
  stop("Paired has to be no|yes")
}

if(is.na(threads)) {
  stop("Threads has to be an integer number")
}

stranded <- if(stranded == "no") 0 else if(stranded == "yes") 1 else 2
paired <- if(paired == "no") 0 else if(paired == "yes") 1 else 2


########################################
## Run dupRadar (have to write own function because dupRadar doesn't let you edit settings for GFF file input etc.)

count <- function(mh, dup) {
    Rsubread::featureCounts(files = bam,
	annot.ext = gtf,
	isGTFAnnotationFile = TRUE,
	GTF.featureType = "CDS",
	GTF.attrType = "ID",
	nthreads = threads,
	isPairedEnd = paired,
	strandSpecific = stranded,
	ignoreDup = dup,
	countMultiMappingReads = mh,)
    }

silencer <- capture.output(counts <- list(mhdup = count(mh = TRUE, dup = FALSE), mhnodup = count(mh = TRUE, dup = TRUE), nomhdup = count(mh = FALSE, dup = FALSE), nomhnodup = count(mh = FALSE, dup = TRUE)))

x <- lapply(counts, function(x) {
    N <- sum(x$stat[, 2]) - x$stat[x$stat$Status == "Unassigned_Unmapped", 2]
    x <- data.frame(gene = rownames(x$counts), width = x$annotation$Length[match(rownames(x$counts), x$annotation$GeneID)], counts = x$counts[, 1], RPK = 0, RPKM = 0)
    x$RPK <- x$counts * (10^3/x$width)
    x$RPKM <- x$RPK * (10^6/N)
    return(x)
    })

x <- data.frame(ID = x[[1]]$gene, geneLength = x[[1]]$width,
    allCountsMulti = x[[1]]$counts, filteredCountsMulti = x[[2]]$counts,
    dupRateMulti = (x[[1]]$counts - x[[2]]$counts)/x[[1]]$counts,
    dupsPerIdMulti = x[[1]]$counts - x[[2]]$counts, RPKMulti = x[[1]]$RPK,
    RPKMMulti = x[[1]]$RPKM, allCounts = x[[3]]$counts, filteredCounts = x[[4]]$counts,
    dupRate = (x[[3]]$counts - x[[4]]$counts)/x[[3]]$counts,
    dupsPerId = x[[3]]$counts - x[[4]]$counts, RPK = x[[3]]$RPK,
    RPKM = x[[3]]$RPKM)


########################################
## Create plots

## duprate vs. expression smooth scatter
png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_drescatter.png"),
    width=1000, height=1000)
duprateExpDensPlot(x, main=basename(bam))
dev.off()

## expression histogram
png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_ehist.png"),
    width=1000, height=1000)
expressionHist(x)
dev.off()

## duprate vs. expression boxplot
png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_drebp.png"),
    width=1000, height=1000)
par(mar=c(10,4,4,2)+.1)
duprateExpBoxplot(x, main=basename(bam))
dev.off()


# 2d density scatter plot
# relating the normalized number of counts per gene (RPK, as a quantification of the gene expression) and the fraction represented by duplicated reads.
# For lowly expressed genes we expect read duplication to happen rarely by chance, while for highly expressed genes - depending on the total sequencing depth - we expect read duplication to happen often.

#par(mfrow=c(1,2))
#duprateExpDensPlot(DupMat=x)

fit <- duprateExpFit(DupMat=x)
cat("duprate at low read counts: ",fit$intercept,"\n",
    "progression of the duplication rate: ",fit$slope,fill=TRUE)
