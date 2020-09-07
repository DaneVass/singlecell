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
num.snps <- nrow(vcf@fix)

# transpose to have cells as rows and SNPs as cols
dat <- t(as.data.frame(vcf@gt))
df <- data.frame(genotype = dat[-1,1:num.snps])

# find rows that 
not.empty <- which(test.row != ".:.:.:.:.:.")
for (i in 1:num.snps){
  geno <- paste("genotype", as.character(i), sep = '.')
  df[,geno] <- sapply(strsplit(as.character(df[,geno]),"/"), `[`, 1)
}

# collapse calls into final annotation
df$final <- as.character(gsub(".:.:.:.:.:.","",apply(df, 1, paste, collapse="")))

# resolve multi SNP rows
table(df$final)
multi.call.rows <- which(nchar(df$final) > 1)
replace <- lapply(multi.call.rows, function(x){
  call <- ifelse(length(unique(strsplit(df$final[x], "")[[1]])) == 1, strsplit(df$final[x], "")[[1]][1], "Unknown")
  df$final[x] <- call
})

# if there is only one unique call label it as such else label it as Unknown
df$final[multi.call.rows] <- unlist(replace)

# rename unknown alleles
rows <- which(df$final == "")
df$final[rows] <- "Unknown"

# rename CD45.2 alleles
rows <- which(df$final == "0")
df$final[rows] <- "CD45.2"

# rename CD45.1 alleles
rows <- which(df$final == "1")
df$final[rows] <- "CD45.1"

ref <- which(df$final == "CD45.2")
alt <- which(df$final == "CD45.1")
unknown <- which(df$final == "Unknown")

# output parsed data to csv
out.df <- df[,"final", drop = F]
write.csv(out.df, outfile)

# report final metrics
message("Filter cellSNP Finished")
message("-------------------------")
message(paste("Processing:", samplename))
message(paste("Output file:", outfile))
message(paste("Number of cells with CD45.2 allele:", length(ref)))
message(paste("Number of cells with CD45.1 allele:", length(alt)))
message(paste("Number of cells with Unknown allele:", length(unknown)))





