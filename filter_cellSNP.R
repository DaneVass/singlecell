#!/usr/bin/env Rscript

# Dane Vassiliadis
# Peter MacCallum Cancer Center
# 7-9-20

# DESCRIPTION
# script to filter the output of cellSNP for ptprc SNPs
# and output a csv compatible with Seurat objects

##------------------------ Setup Environment ------------------------##
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("stringr"))

##------------------------ Get Options ------------------------##
option_list <- list(
  # Basic options
  make_option(c("-i", "--infile"), type="character", dest="infile", default = NULL,
              help="REQUIRED: Path to the input vcf file outputted by cellSNP", metavar="input file"),
  make_option(c("-o", "--outdir"), type="character", dest="outdir", default="./",
              help="Folder to output analysis files."),
  make_option(c("-n", "--samplename"), type="character", dest="samplename", default="seurat_pipeline_default",
              help="Name of sample. Default = seurat_pipeline_default"))

#-----------------------
# Parse options
#-----------------------
options <- parse_args(OptionParser(option_list=option_list))

#-----------------------
# Run pipeline
#-----------------------
if(is.null(options$infile)){
  message("no input file given. quitting")
  quit(status = 1)
} else {
  infile = options$infile
}
outdir = options$outdir
samplename = options$samplename
outfile = file.path(outdir, paste(samplename,"_Ptprc_genotype.csv", sep = ''))


message("Filter cellSNP Start")
message("-------------------------")
message(paste("Processing:", samplename))
message(paste("Input file:", basename(infile)))
message(paste("Output file:", outfile))


# Filter cellSNP  output
vcf <- vcfR::read.vcfR(infile, verbose = T)

# transpose to have cells as rows and SNPs as cols
dat <- t(as.data.frame(vcf@gt))
df <- data.frame(cell = rownames(dat), genotype = dat[,1])

# how many cells do not have genotype assigned
remove <- c(1,which(df$genotype == ".:.:.:.:.:."))
df <- df[-remove,,drop = F]
df$background <- sapply(strsplit(as.character(df$genotype),":"), `[`, 1)
ref <- which(df$background == "0/0")
alt <- which(df$background == "1/1")

# output parsed data to csv
write.csv(df, outfile)

# report final metrics
message("Filter cellSNP Finished")
message("-------------------------")
message(paste("Processing:", samplename))
message(paste("Output file:", outfile))
message(paste("Number of cells with REF allele:", length(ref)))
message(paste("Number of cells with ALT allele:", length(alt)))
message(paste("Number of cells with NO assigned allele:", length(remove)))





