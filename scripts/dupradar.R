#!/usr/bin/env Rscript

########################################
##
## dupRadar shell script
## call dupRadar R package from the shell for
## easy integration into pipelines
##
## Holger Klein & Sergi Sayols
##
## https://sourceforge.net/projects/dupradar/
##
## input:
## - _duplicate marked_ bam file
## - gtf file
## - parameters for duplication counting routine:
##   stranded, paired, outdir, threads.
##
########################################

library(dupRadar)

####################
##
## get name patterns from command line
##
args   <- commandArgs(TRUE)

## the bam file to analyse
bam <- args[1]
## usually, same GTF file as used in htseq-count
gtf <- gsub("gtf=","",args[2])
## no|yes|reverse
stranded <- gsub("stranded=","",args[3])
## is a paired end experiment
paired   <- gsub("paired=","",args[4])
## output directory
outdir   <- gsub("outdir=","",args[5])
## number of threads to be used
threads  <- as.integer(gsub("threads=","",args[6]))

if(length(args) != 6) {
  stop (paste0("Usage: ./dupRadar.R <file.bam> <genes.gtf> ",
               "<stranded=[no|yes|reverse]> paired=[yes|no] ",
               "outdir=./ threads=1"))
}

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

## end command line parsing
##
########################################

########################################
##
## analyze duprates and create plots
##
cat("Processing file ", bam, " with GTF ", gtf, "\n")

## calculate duplication rate matrix
dm <- analyzeDuprates(bam,
                      gtf,
                      stranded,
                      (paired == "yes"),
                      threads)

## produce plots

## duprate vs. expression smooth scatter
png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_drescatter.png"),
    width=1000, height=1000)
duprateExpDensPlot(dm, main=basename(bam))
dev.off()

## expression histogram
png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_ehist.png"),
    width=1000, height=1000)
expressionHist(dm)
dev.off()

## duprate vs. expression boxplot
png(file=paste0(outdir,"/",gsub("(.*)\\.[^.]+","\\1",basename(bam)),"_dupRadar_drebp.png"),
    width=1000, height=1000)
par(mar=c(10,4,4,2)+.1)
duprateExpBoxplot(dm, main=basename(bam))
dev.off()
