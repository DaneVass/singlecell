#!/usr/bin/env/Rscript

# seurat_integration.R
# Dane Vassiliadis
# Peter MacCallum Cancer Center
# 23-6-20

# DESCRIPTION
# Shell R Pipeline to automate single cell sample integration using Seurat CCA analysis

##------------------------ Setup Environment ------------------------##
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("sctransform"))
suppressPackageStartupMessages(library("DropletUtils"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("stringr"))

##------------------------ Get Options ------------------------##
option_list <- list(
    
    # Basic options
    make_option(c("-i", "--inputdir"), type="character", dest="inputdir", default = NULL,
            help="REQUIRED: Path to the input directory containing Seurat objects for integration. It is easiest to symlink all Seurat objects that you want to integrate into this folder", metavar="path"),
    make_option(c("-o", "--outdir"), type="character", dest="outdir", default="./",
            help="Folder to output integrated Seurat object. Will be created if it does not exist"),
    make_option(c("-n", "--samplename"), type="character", dest="samplename", default="seurat_integration_CCA_default",
            help="Name of sample. [default]"),
    make_option(c("--meta"), type="character", dest="meta",
            help="Metadata to add to the sample. Usually a table of cellID with metadata columns. Must contain same cell IDs as the sample."),
    make_option(c("--norm"), type="character", dest="norm", default="sctransform",
            help="Normalisation method to use [default]"),

    # Integration options
    make_option(c("-f", "--features"), type="integer", default=3000, dest="features",
            help="Number of Integration features to select. [default]")
    make_option(c("--norm"), type="character", dest="norm", default="sctransform",
            help="Normalisation method to use [default]"),
    
    # Other options
    make_option(c("-d", "--dims"), type="integer", default=50, dest="dims",
            help="Upper limit of PCA dimensions to use for clustering and t-SNE/UMAP projection. Default = 30"),
    make_option(c("--regressvars"), type="character", dest="regressvars", default='c("percent.mito", "S.Score", "G2M.Score")', metavar = "vars to regress",
            help="Variables to regress out during sample normalisation. Default = 'percent.mito, S.Score, G2M.Score'"),

    # Runtime options
    make_option("--multithreading", type="store_true", default=TRUE, dest="multithreading", 
            help="Use multithreading via the future package to speedup computation"),
    make_option(c("-t", "--threads"), type="integer", default=8, dest="threads", 
            help="Desired number of CPU threads to use [default]"),
    make_option(c("-m", "--memory"), type="integer", default=16, dest="memory", 
            help="Desired RAM in GB. E.g. for 16GB: -m 16 [default]"),
    make_option(c("--seed"), type="integer", default=10101, dest="seed",
            help="Random seed to use [default]")        

)

#-----------------------
# Parse options
#-----------------------
options <- parse_args(OptionParser(option_list=option_list))

# check inputs
if(is.null(options$inputdir)){
    message("no input directory given. quitting")
    quit(status = 1)
} else {
    input.dir = file.path(options$inputdir)
}

# collect input files
input.files <- list.files(input.dir)




#-----------------------
# Run integration pipeline
#-----------------------

message("---------------------------")
message("Seurat integration pipeline")
message("---------------------------")
message(paste("Input directory:", counts.dir))
message(paste("Found the following files:", input.files))
message(paste("Saving integrated object to:", options$outdir))
message(paste("Output filename:", options$samplename))
message(paste("Integration method:", "CCA")) # TODO: implement other integration methods such as harmony and LIGER

message(paste("Variables to regress during normalisation:", toString(regressVars)))
message(paste("PCA dimensions to use:", toString(dims)))
message(paste("Clustering resolutions to use:", toString(resolution)))
message(paste("Number of integration features to select:", options$features))
message(paste("Random seed:", seed.use))

# setup output dirs
output.dir <- file.path(outdir)
dir.create(output.dir, showWarnings = F)
  
message("")
message("Reading in Seurat datasets")

obj.list <- lapply(input.files, function(x){
    readRDS(file.path(x))
    })
str(obj.list)

# Set max file size options 
mem = options$memory * 1000
options(future.globals.maxSize = mem * 1024^2)

if (isTRUE(options$multithreading)){
    require(future)
    # change the current plan to access parallelization
    plan("multiprocess", workers = options$threads)
    plan()
}

# make sure each sample is named for later
for i in 1:length(obj.list){
    obj.list[i]$samplename <- as.character(input.files[i]) 
}

# reset default assay
for i in 1:length(obj.list){
    DefaultAssay(obj.list[i]) <- "SCT"
}

# CCA integration workflow
integration.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = options$features)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = integration.features, verbose = T)
integration.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = integration.features, verbose = T)
obj.integrated <- IntegrateData(anchorset = integration.anchors, normalization.method = "SCT", verbose = T)

obj.integrated <- RunPCA(object = obj.integrated, verbose = T, seed.use = options$seed, npcs = 100)
ElbowPlot(obj.integrated, ndims = 100)
ggsave(paste(options$samplename, "integration.pdf", sep = "_")
obj.integrated <- RunUMAP(object = obj.integrated, dims = 1:options$dims)

# save output file
saveRDS(obj.integrated, file = file.path(paste(options$samplename, "integrated.rds", sep = "_")))