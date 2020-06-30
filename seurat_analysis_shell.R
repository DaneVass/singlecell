#!/usr/bin/env Rscript

# Dane Vassiliadis
# Peter MacCallum Cancer Center
# 10-6-20

# DESCRIPTION
# Pipeline to automate basic 10X sample processing and visualisation using
# Seurat 

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
  make_option("--object", type="character", dest="object", default = NULL, 
      help="OPTIONAL: Path to the input Seurat RDS object.", metavar="path"),
  make_option(c("-o", "--outdir"), type="character", dest="outdir", default="./",
      help="Folder to output analysis files. Will be created if it does not exist"),
  make_option(c("-n", "--samplename"), type="character", dest="samplename", default="seurat_pipeline_default",
      help="Name of sample. Default = seurat_pipeline_default"),
  make_option(c("-s", "--species"), type="character", dest="species", default="human",
      help="Species of sample. Can be human or mouse. Default = human"),
  make_option(c("--project"), type="character", dest="project", default="scRNAseq",
      help="Project name sample belongs to. Default = scRNAseq"),
  make_option(c("--meta"), type="character", dest="meta",
      help="Metadata to add to the sample. Usually a table of cellID with metadata columns. Must contain same cell IDs as the sample."),
  make_option(c("--norm"), type="character", dest="norm", default="sctransform",
      help="Normalisation method to use [default]"),

  # Filtering options
  make_option(c("-m", "--mito"), type="integer", default=10, dest="mitocutoff",
    help="Threshold for filtering droplets based on mitochondrial read percentage. Default = 10"),
  make_option("--featuresup", type="integer", default=10000, dest="featuresup",
    help="Upper threshold for filtering droplets based on total detected features. Default = 10000"),
  make_option("--featuresdn", type="integer", default=100, dest="featuresdn",
    help="Lower threshold for filtering droplets based on total detected features. Default = 100"),
  make_option("--countsup", type="integer", default=100000, dest="countsup",
    help="Upper threshold for filtering droplets based on total read counts. Default = 100000"),
  make_option("--countsdn", type="integer", default=1000, dest="countsdn",
    help="Lower threshold for filtering droplets based on total detected features. Default = 1000"),
    
  # Other options
  make_option(c("-d", "--dims"), type="integer", default=30, dest="dims",
    help="Upper limit of PCA dimensions to use for clustering and t-SNE/UMAP projection. Default = 30"),
  make_option(c("--regressvars"), type="character", dest="regressvars", default='c("percent.mito", "S.Score", "G2M.Score")', metavar = "vars to regress",
    help="Variables to regress out during sample normalisation. Default = 'percent.mito, S.Score, G2M.Score'"),
    
  # Runtime options
  make_option(c("-t", "--threads"), type="integer", default=4, dest="threads", 
    help="Desired number of CPU threads to use [default]"),
  make_option(c("--seed"), type="integer", default=10101, dest="seed",
    help="Random seed to use [default]")        
)

#-----------------------
# Parse options
#-----------------------
options <- parse_args(OptionParser(option_list=option_list))

#-----------------------
# Run pipeline
#-----------------------
if(is.null(options$countsdir)){
  if(is.null(options$object)){
    message("no input directory or object given. quitting")
    quit(status = 1)
  } else {
      obj = options$object
      counts.dir = NULL
} else {
  counts.dir = options$countsdir
  obj = NULL
}
outdir = options$outdir 
norm = options$norm
mitoCutoff = options$mitocutoff
featuresLower = options$featuresdn
featuresUpper = options$featuresup
countsLower = options$countsdn
countsUpper = options$countsup 
dims = 1:options$dims 
verbose = options$verbose
resolution <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
samplename = options$samplename
regressVars = options$regressvars 
species = options$species 
project = options$project
meta = options$meta 
seed.use = options$seed

message("-------------------------")
message("Seurat scRNAseq pipeline")
message("-------------------------")
message(paste("Processing:", samplename))
if(!is.null(options$countsdir)){message(paste("Input directory:", counts.dir))}
if(!is.null(options$object)){message(paste("Input object:", obj))}
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
if(!is.null(options$countsdir)){
  message("Importing count data")
  counts <- Read10X(counts.dir, strip.suffix = T)
  sample <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = project, meta.data = meta)
}

# import Seurat object
if(is.null(options$obj)){
    sample <- readRDS(file = obj)
}

#-----------------------
# Sample QC
#-----------------------
message("Performing sample QC")
mito.pattern <- "^mt-"
if(species == "human" | species == "Human" | species == "Hs" | species == "hs"){
  mito.pattern <- "^MT-"
}

# calculate mito percentage
message("Calculating mito %")
mito.genes <- grep(pattern = mito.pattern, x = rownames(x = sample@assays$RNA), value = TRUE)
percent.mito <- Matrix::colSums(sample@assays$RNA[mito.genes, ])/Matrix::colSums(sample@assays$RNA)*100
sample <- AddMetaData(object = sample, metadata = percent.mito, col.name = "percent.mito")

# plot QC metrics
message("Plotting QC metrics - Pre-filter")
pdf(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_pre-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T)
dev.off()

# Feature Scatter plots
pdf(file.path(qc.dir,paste(samplename, "_count_mito_QCscatter-pre-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

pdf(file.path(qc.dir,paste(samplename, "_count_feature_QCscatter-pre-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# subset based on mito count and feature cutoffs
message("Filtering low quality cells")
message("Starting dimensions")
message(print(dim(sample)))
sample <- subset(sample, cells = colnames(sample)[which(sample$percent.mito < mitoCutoff)])
sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA > countsLower)])
sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA < countsUpper)])
sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA > featuresLower)])
sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA < featuresUpper)])

message("")
message("Filtered dimensions")
message(print(dim(sample)))

# plot QC metrics
message("Plotting QC metrics - Post-filter")
pdf(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_post-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T)
dev.off()

# Feature Scatter plots
pdf(file.path(qc.dir,paste(samplename, "_count_mito_QCscatter-postfilter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()
pdf(file.path(qc.dir,paste(samplename, "_count_feature_QCscatter-post-filter.pdf", sep = '')), width = 9, height = 7, useDingbats = F)
FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


#-----------------------
# Normalisation #1
#-----------------------
message(paste0("Normalising counts using: ", norm))
# SCTransform workflow
message("Regressing following variables: percent.mito")

if(norm == "sctransform") {
  sample <- SCTransform(sample, vars.to.regress = "percent.mito", verbose = T, conserve.memory = F, return.only.var.genes = F, seed.use = seed.use)
}

# log normalisation workflow
if(norm == "lognorm") {
  sample <- NormalizeData(sample, normalization.method = "LogNormalize",
                            vars.to.regress = regress.vars,
                            verbose = F, conserve.memory = F)
}

#----------------------------------------------
# Seurat workflow without cell cycle regression
#----------------------------------------------
message("Seurat workflow without cell cycle regression")

# run PCA
message("Running PCA")
sample <- RunPCA(sample, npcs = 100, verbose = T, seed.use = seed.use)

# Run tSNE & UMAP
message("Running UMAP")
sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap",tail(dims,1), sep='_'))

message("Running TSNE")
sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne",tail(dims,1), sep='_'))

message("Finding nearest neighbors")
sample <- FindNeighbors(object = sample, verbose = T, dims = dims)

message("Running cluster sequence")
# run a sequence of cluster resolutions
for (i in 1:length(resolution)){
  sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
}

message("Saving normalised dataset pre cell cycle regression")
# save normalised and clustered dataset
saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_normalised_no_cc_regression.rds", sep = "")))

#----------------------------------------------
# Seurat workflow with cell cycle regression
#----------------------------------------------

#-----------------------
# Cell cycle scoring
#-----------------------

message("Performing Cell cycle scoring")
# perform cell cycle analysis
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

if (species == "mouse"){
  message("Species is mouse. Converting human to mouse gene symbols")
  convertHumanGeneList <- function(x){
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
      values = x , mart = human, attributesL = c("mgi_symbol"),
      martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2]) 
    return(humanx)
  }
  s.genes <- convertHumanGeneList(s.genes)
  g2m.genes <- convertHumanGeneList(g2m.genes)
}

sample <- CellCycleScoring(sample, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, search = F)

# add in cell cycle difference scoring
sample$CC.Difference <- sample$S.Score - sample$G2M.Score

# run PCA on cell cycle components
sample <- RunPCA(sample, features = c(s.genes, g2m.genes))
p <- DimPlot(sample, reduction = "pca", group.by = 'Phase', label = T, label.size = 5)
ggsave(file.path(plot.dir,paste(samplename, "_PCA_cell_cycle_anno.pdf", sep = '')), width = 9, height = 7)

#-----------------------
# Normalisation #2
#-----------------------
message("Normalisation #2")
message(paste("Normalising counts using:", norm))
# SCTransform workflow
message("Regressing following variables:")

regressVars <- str_split(regressVars, ',')[[1]]
print(regressVars)

if(norm == "sctransform") {
  sample <- SCTransform(sample, vars.to.regress = regressVars, verbose = T, conserve.memory = F, return.only.var.genes = F)
}

# log normalisation workflow
if(norm == "lognorm") {
  sample <- NormalizeData(sample, normalization.method = "LogNormalize",
                            vars.to.regress = regress.vars,
                            verbose = F, conserve.memory = F)
}

#--- Dimensional reduction & clustering
message("Running dimensional reduction algos")

# always calculate 100 PCs
sample <- RunPCA(object = sample, verbose = T, npcs = 100, seed.use = seed.use)
message("Plotting...")

# Elbow plot
percent.var <- Stdev(sample)
elbow.dim <- PCAtools::findElbowPoint(percent.var)
message(paste("Predicted elbow PC =", elbow.dim))
ElbowPlot(sample, ndims = 50) +
  geom_vline(aes(xintercept = elbow.dim), color = "red") +
  labs(title = paste("Elbow plot -",samplename), subtitle = paste("Predicted elbow PC =", elbow.dim))
ggsave(file.path(plot.dir,paste(samplename,"_elbow_plot.pdf", sep = '')), width = 8, height = 6)

# Run tSNE & UMAP on 1 to elbow dimensions
message("Running UMAP")
sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name="umap_elbow")
message("Running TSNE")
sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name="tsne_elbow")

# find neighbors based on 1 to elbow dimensions
sample <- FindNeighbors(object = sample, verbose = T, dims = 1:elbow.dim)

# Run tSNE & UMAP on 1 to dims
message("Running UMAP")
sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap",tail(dims,1), sep='_'))
message("Running TSNE")
sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne",tail(dims,1), sep='_'))

# run a sequence of cluster resolutions for elbow dim
for (i in 1:length(resolution)){
  sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
}

# save normalised and clustered dataset
message("Saving normalised dataset post cell cycle regression")
saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_cc_regressed_normalised.rds", sep = "")))

message("performing final plots on normalised dataset")

# PCA plot reduction on clusters
p <- DimPlot(object = sample, reduction = "pca", pt.size = 0.1, label = T, label.size = 5, dims = c(head(dims,1),tail(dims,1))) +
  NoLegend() +
  ggtitle(paste('PCA -', samplename))
p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
ggsave(file.path(plot.dir,paste(samplename, "_PCA.pdf", sep = '')), width = 9, height = 7)

# UMAP reduction on clusters - elbow
p <- DimPlot(object = sample, reduction = "umap_elbow", pt.size = 0.1, label = T, label.size = 5) +
  NoLegend() +
  ggtitle(paste('UMAP -', samplename))
p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
ggsave(file.path(plot.dir,paste(samplename, "_UMAP_elbow.pdf", sep = '')), width = 9, height = 7)

# t-SNE reduction on clusters - elbow
p <- DimPlot(object = sample, reduction = "tsne_elbow", pt.size = 0.1, label = T, label.size = 5) +
  NoLegend() +
  ggtitle(paste('t-SNE -', samplename))
p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
ggsave(file.path(plot.dir,paste(samplename, "_t-SNE_elbow.pdf", sep = '')), width = 9, height = 7)

# # UMAP reduction on clusters - dims
# p <- DimPlot(object = sample, reduction = paste("umap",tail(dims, 1)), sep='_'), pt.size = 0.1, label = T, label.size = 5) +
#   NoLegend() +
#   ggtitle(paste('UMAP -', samplename))
# p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
# ggsave(file.path(plot.dir,paste(samplename, paste("_UMAP_",tail(dims, 1)), sep=''), ".pdf", sep = '')), width = 9, height = 7)
  
# # t-SNE reduction on clusters - dims
# p <- DimPlot(object = sample, reduction = paste("tsne",dims), sep='_'), pt.size = 0.1, label = T, label.size = 5) +
#   NoLegend() +
#   ggtitle(paste('t-SNE -', samplename))
# p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
# ggsave(file.path(plot.dir,paste(samplename, paste("_t-SNE_",dims), sep=''), ".pdf", sep = '')), width = 9, height = 7)
  
message(paste("Finished", samplename, "analysis. Enjoy your data!"))

sessionInfo()
