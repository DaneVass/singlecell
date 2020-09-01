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
    
)

#-----------------------
# Parse options
#-----------------------
options <- parse_args(OptionParser(option_list=option_list))

#-----------------------
# Run integration pipeline
#-----------------------
if(is.null(options$inputdir)){
    message("no input directory given. quitting")
    quit(status = 1)
} else {
    input.dir = file.path(options$inputdir)
}

message("-------------------------")
message("Seurat scRNAseq pipeline")
message("-------------------------")
message(paste("Processing:", samplename))
message(paste("Input directory:", counts.dir))
message(paste("Output directory:", outdir))
message(paste("Species:", species))
message(paste("Project:", project))
message(paste("Mitochondrial % cutoff:", mitoCutoff))
message(paste("# features cutoff:", featuresLower, featuresUpper))
message(paste("# read counts cutoff:", countsLower, countsUpper))
message(paste("Variables to regress during normalisation:", toString(regressVars)))
message(paste("PCA dimensions to use:", toString(dims)))
message(paste("Clustering resolutions to use:", toString(resolution)))
#message(paste("Metadata:", meta))
message(paste("Random seed:", seed.use))

  # setup output dirs
  output.dir <- file.path(outdir)
  dir.create(output.dir, showWarnings = F)

  plot.dir <- file.path(outdir, "plots")
  dir.create(plot.dir, showWarnings = F)

  qc.dir <- file.path(outdir, "qc")
  dir.create(qc.dir, showWarnings = F)

  Seurat.obj.dir <- file.path(outdir, "seurat_obj")
  dir.create(Seurat.obj.dir, showWarnings = F)

  message("")
  message("Importing count data")
  counts <- Read10X(counts.dir, strip.suffix = T)
  sample <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = project, meta.data = meta)



# Set max file size options 
options(future.globals.maxSize = 16000 * 1024^2)

# make sure each sample is named for later
KRAS.T0$samplename <- "KRAS.T0"
KRAS.14$samplename <- "KRAS.14"
KRAS.17$samplename <- "KRAS.17"

# reset default assay
DefaultAssay(KRAS.T0) <- "SCT"
DefaultAssay(KRAS.14) <- "SCT"
DefaultAssay(KRAS.17) <- "SCT"

# CCA integration workflow
KRAS.list <- list(KRAS.T0, KRAS.14, KRAS.17) 
KRAS.features <- SelectIntegrationFeatures(object.list = KRAS.list, nfeatures = 6000)
KRAS.list <- PrepSCTIntegration(object.list = KRAS.list, anchor.features = KRAS.features, verbose = T)
KRAS.anchors <- FindIntegrationAnchors(object.list = KRAS.list, normalization.method = "SCT", anchor.features = KRAS.features, verbose = T)
KRAS.integrated <- IntegrateData(anchorset = KRAS.anchors, normalization.method = "SCT", verbose = T)
KRAS.integrated@meta.data
KRAS.integrated <- RunPCA(object = KRAS.integrated, verbose = T)
KRAS.integrated <- RunUMAP(object = KRAS.integrated, dims = 1:50)


## CCA integration
```{r CCA integration}
# We will use the seurat approach to integrate the MLL_T0 datasets and perform batch correction
# https://satijalab.org/seurat/v3.1/integration.html

options(future.globals.maxSize = 32000 * 1024^2)

DefaultAssay(P1) <- "SCT"
DefaultAssay(P2) <- "SCT"

# relabel samples based on hashtag information
P1$samplename <- P1$HTO_classification
P2$samplename <- P2$HTO_classification

# CCA integration workflow
TRAM.list <- list(P1, P2) 
TRAM.features <- SelectIntegrationFeatures(object.list = TRAM.list, nfeatures = 3000)
TRAM.list <- PrepSCTIntegration(object.list = TRAM.list, anchor.features = TRAM.features, verbose = T)
TRAM.anchors <- FindIntegrationAnchors(object.list = TRAM.list, normalization.method = "SCT", anchor.features = TRAM.features, verbose = T)
TRAM.int <- IntegrateData(anchorset = TRAM.anchors, normalization.method = "SCT", verbose = T)
TRAM.int@meta.data
TRAM.int <- RunPCA(object = TRAM.int, verbose = T, seed.use = 10101, npcs = 100)
ElbowPlot(TRAM.int, ndims = 100)
TRAM.int <- quickUMAPrerun(TRAM.int)
DimPlot(TRAM.int)
DimPlot(TRAM.int, group.by = "samplename")
saveRDS(TRAM.int, file = file.path(TRAM.results.dir, "P1_P2_integrated.rds"))