#!/usr/bin/env/Rscript

# seurat_integration.R
# Dane Vassiliadis
# Peter MacCallum Cancer Center
# 23-6-20

# DESCRIPTION
# Pipeline to automate basic 10X sample processing and visualisation using
# Seurat CCA analysis

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
    make_option(c("-i", "--countsdir"), type="character", dest="countsdir", default = NULL,
        help="REQUIRED: Path to the input directory containing counts matrix, barcodes, and features files", metavar="path"),
 25   make_option(c("-o", "--outdir"), type="character", dest="outdir", default="./",
 26       help="Folder to output analysis files. Will be created if it does not exist"),
 27   make_option(c("-n", "--samplename"), type="character", dest="samplename", default="seurat_pipeline_default",
 28       help="Name of sample. Default = seurat_pipeline_default"),
 29   make_option(c("-s", "--species"), type="character", dest="species", default="human",
 30       help="Species of sample. Can be human or mouse. Default = human"),
 31   make_option(c("--project"), type="character", dest="project", default="scRNAseq",
 32       help="Project name sample belongs to. Default = scRNAseq"),
 33   make_option(c("--meta"), type="character", dest="meta",
 34       help="Metadata to add to the sample. Usually a table of cellID with metadata columns. Must contain same cell IDs as the sa    mple."),
 35   make_option(c("--norm"), type="character", dest="norm", default="sctransform",
 36       help="Normalisation method to use [default]"),
 37
 38   # Filtering options
 39   make_option(c("-m", "--mito"), type="integer", default=10, dest="mitocutoff",
 40     help="Threshold for filtering droplets based on mitochondrial read percentage. Default = 10"),
 41   make_option("--featuresup", type="integer", default=10000, dest="featuresup",
 42     help="Upper threshold for filtering droplets based on total detected features. Default = 10000"),
 43   make_option("--featuresdn", type="integer", default=100, dest="featuresdn",
 44     help="Lower threshold for filtering droplets based on total detected features. Default = 100"),
 45   make_option("--countsup", type="integer", default=100000, dest="countsup",
 46     help="Upper threshold for filtering droplets based on total read counts. Default = 100000"),
 47   make_option("--countsdn", type="integer", default=1000, dest="countsdn",
 48     help="Lower threshold for filtering droplets based on total detected features. Default = 1000"),


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