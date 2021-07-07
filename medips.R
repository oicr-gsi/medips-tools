library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MEDIPSData)
library(optparse)

# command line options
option_list = list(
  make_option(c("-d", "--basedir"), type="character", default=NULL, help="base directory", metavar="character"),
  make_option(c("-i", "--bamfile"), type="character", default=NULL, help="input bam file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="output directory to be created", metavar="character"),
  make_option(c("-g", "--BSgenome"), type="character", default="BSgenome.Hsapiens.UCSC.hg19", help="genome", metavar="character"),
  make_option(c("-u", "--uniq"), type="numeric", default=1e-3, help="uniq", metavar="numeric"),
  make_option(c("-e", "--extend"), type="numeric", default=300, help="extend", metavar="numeric"),
  make_option(c("-s", "--shift"), type="numeric", default=0, help="shift", metavar="numeric"),
  make_option(c("-w", "--ws"), type="numeric", default=100, help="ws", metavar="numeric"),
  make_option(c("-n", "--samplename"), type="character", default=NULL, help="ws", metavar="character"),
  make_option(c("-f", "--cigarFlag"), action="store_true", default=FALSE, help="simpleCigarFlag"), # add a flag enabling or disabling simple cigar flag setting
  make_option(c("-x", "--covX"), type="numeric", default=1, help="minimum reads supporting CpGs", metavar="numeric") # min cov
  
)

# get options
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# set better variable names
basedir <- opt$basedir
bamfile <- opt$bamfile
outdir <- opt$outdir
BSgenome <- opt$BSgenome
uniq <- opt$uniq
extend <- opt$extend
shift <- opt$shift
ws <- opt$ws
samplename <- opt$samplename
simpleCigarFlag <- opt$cigarFlag
covX <- opt$covX


# set the working dir
setwd(basedir)

# set chromosomes
chr.select=paste0("chr",c(1:22,"X","Y"))

# set sample name if it exists
if (exists("samplename")) {
	sample_name <- samplename
} else {
	sample_name <- basename(bamfile)
}

# print options to output
print("Running MEDIPS with the following options:")
print(opt)

# create output if does not exist
if(!dir.exists(outdir)){
  dir.create(outdir)
}

# run MEDIPS and extract counts
MeDIP.set=MEDIPS.createSet(file=bamfile, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=chr.select, simpleCigar = simpleCigarFlag)

# coupling
CS <- MEDIPS.couplingVector(pattern="CG", refObj=MeDIP.set)

# get saturation metrics
sr <- MEDIPS.saturation(file=bamfile, BSgenome=BSgenome, uniq=uniq, extend=extend, shift=shift, window_size=ws, chr.select=chr.select, simpleCigar = simpleCigarFlag)

# Coverage
cr <- MEDIPS.seqCoverage(file=bamfile, pattern="CG", BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, chr.select=chr.select, simpleCigar = simpleCigarFlag)
# compute the % CpG sites covered at Nx coverage
total.CpG <- length(cr$cov.res)
CpG.with.coverage <- length(ch1$cov.res[ch1$cov.res >= covX]) 
frac.CpG.coverage <- CpG.with.coverage/total.CpG

# CpG enrichment
er <- MEDIPS.CpGenrich(file=bamfile, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, chr.select=chr.select)

# write out counts
write.table(data.frame(MeDIP.set@genome_count), paste0(outdir, "/genome_count.txt"), row.names=F, quote=F, col.names=F)

# write out saturation metrics
saturation_df <- data.frame(maxEstCorReads=sr$maxEstCor[1], maxEstCor=sr$maxEstCor[2], maxTruCorReads=sr$maxTruCor[1], maxTruCor=sr$maxTruCor[2])
rownames(saturation_df) <- sample_name
write.table(saturation_df, file=paste0(outdir, "/saturation_metrics.txt"), sep="\t", row.names=T, quote=F, col.names=NA)

# write out coverage metrics
coverage_df <- data.frame(numberReadsCG=cr$numberReads, numberReadsWOCG=cr$numberReadsWO, numberCpGsitesInGenome=total.CpG, numberCpGSitesWithCov=CpG.with.coverage, fractionCpGCov=frac.CpG.coverage)
rownames(coverage_df) <- sample_name
write.table(coverage_df, file=paste0(outdir, "/coverage_counts.txt"), sep="\t", row.names=T, quote=F, col.names=NA)

# write out enrichment metrics
enrichment_df <- t(data.frame(unlist(er)))
rownames(enrichment_df) <- sample_name
write.table(enrichment_df, file=paste0(outdir, "/enrichment_data.txt"), sep="\t", row.names=T, quote=F, col.names=NA)

# get windows
no.window <- NULL
for(i in 1:length(chr.select)){
  no.window <- c(no.window, ceiling(MeDIP.set@chr_lengths[i]/ws))
}
window.last <- NULL
for(i in 1:length(chr.select)){
  window.last <- c(window.last, sum(no.window[1:i]))
}
window.first <- window.last-no.window+1

# write out window co-ordinates
write.csv(data.frame(chr=chr.select, window.first, window.last, no.window=no.window), paste0(outdir,"/MEDIPS_window_per_chr.csv"), row.names=F)

# write wig
MEDIPS.exportWIG(Set=MeDIP.set, file=paste0(outdir,"/medips.wig"), format="rpkm", descr="")

# write saturation metrics
#png(file=paste0(outdir, "/saturation.png"), units="px", height=1000, width=1000)
 #MEDIPS.plotSaturation(sr)
#dev.off()

# write coverage metrics
#png(file=paste0(outdir, "/coverage.png"), units="px", height=1000, width=1000)
 #MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="hist", t=15, main="Sequence pattern coverage, histogram")
#dev.off()

# write calibration plot
#png(file=paste0(outdir, "/calibration.png"), units="px", height=1000, width=1000)
# MEDIPS.plotCalibrationPlot(CSet=CS, main="Calibration Plot", MSet=MeDIP.set, rpkm=TRUE, xrange=TRUE)
#dev.off()
