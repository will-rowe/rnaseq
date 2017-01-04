
###
# PER BASE GENOME COVERAGE
###

# samtools flagstat for quick stats on alignment
$ samtools flagstat alignment.sorted.bam

# create a bed file for your region of interest (e.g. the chromosome)
$ printf 'D23580_liv_chro\t1\t4879402\tD23580.chrom\t0\n' > D23580.chrom.bed

# calculate coverage stats with bedtools
$ bedtools coverage -a D23580.chrom.bed -b alignment.sorted.bam -d | gzip > alignment.coverage.tsv.gz

# number of bases in region of interest (CP002487.fasta)
$ BaseNo=`gunzip -c alignment.coverage.tsv.gz | wc -l`

# find number of bases with coverage > 30 and calculate percentage
$ Cov=`gunzip -c alignment.coverage.tsv.gz | awk '$7>30 {print}' | wc -l`
$ bc -l<<<$Cov/$BaseNo*100


###
# PLOT COVERAGE INFO FROM BAM FILE
###

## A. If you have a target/suspect region, you could try something like this:
# calculate coverage stats per base with bedtools
$ bedtools coverage -a suspect_region.bed -b alignment.sorted.bam -d | gzip > alignment.coverage.tsv.gz

# grab Chr, Locus and Depth
$ gunzip -c alignment.coverage.tsv.gz | cut -f 1,6,7 - > alignment.cov_4_R.tsv

# example R script:
library(reshape)
library(lattice, pos=10)
chrom <- read.table("alignment.cov_4_R.tsv", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
chrom <- rename(chrom,c(V1="Chr", V2="Locus", V3="Depth"))
trellis.device(device='png', filename='file.png')
my.plot <- xyplot(Depth ~ Locus, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), data=chrom, main="coverage by locus - region of interest")
print(my.plot)
dev.off()

## B. For whole genomes, you could smooth the results by taking averages from windows of 1k bases
$ bedtools makewindows -b D23580.chrom.bed -w 1000 > D23.1k.bed
$ bedtools coverage -hist -a D23.1k.bed -b alignment.sorted.bam | grep ^all > alignment.hist.all.txt

## C. Another option is to calculate the cumulative distribution describing the fraction of targeted bases that were covered by >x reads
# use the bedtools coverage histogram
$ bedtools coverage -hist -a suspect_region.bed -b alignment.sorted.bam | grep ^all > alignment.hist.txt

# use GNU parallel for multiple samples
$ find ./ -type f -name '*bam' | parallel 'bedtools coverage -hist -a D23.1k.bed -b {} | grep ^all > {}.hist.all.txt'

# use an R script to visualise coverage over regions of interst / whole genome
## this example is by Stephen Turner
## http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html

# Get a list of the bedtools output files you'd like to read in
print(files <- list.files(pattern="all.txt$"))
print(labels <- paste("Sample-", gsub("\\.sorted.bam.hist.all.txt", "", files, perl=TRUE), sep=""))
# Create lists to hold coverage and cumulative coverage for each alignment and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}
# Pick colours
library(RColorBrewer)
cols <- brewer.pal(length(cov), "Dark2")
# Save the graph to a file
png("coverage-plots.png", h=1000, w=1000, pointsize=20)
# Create plot area, but do not plot anything. Add gridlines and axis labels.
plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth", ylab="Fraction of genome \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))
# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=cols[i])
# Add a legend using the nice sample labeles rather than the full filenames.
par(xpd=TRUE)
legend("topright", legend=labels, col=cols, lty=1, lwd=4)
dev.off()
