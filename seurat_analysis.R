# runSeuratPipeline
# Dane Vassiliadis
# 8-6-20

# Pipeline to automate basic 10X sample processing and visualisation using
# Seurat


#-----------------------
# Required libs
#-----------------------
require(Seurat)
require(sctransform)
require(DropletUtils)
require(PCAtools)

#-----------------------
# Main function
#-----------------------
runSeuratPipeline <- function(counts.dir, outdir = "./", norm = "sctransform",
  mitoCutoff = 20,featuresLower = 100L,
  featuresUpper = 10000L,countsLower = 1000L,
  countsUpper = 100000L, dims = 1:30, verbose = T,
  resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
  samplename = "seurat_pipeline_default",
  regressVars = c("percent.mito", "S.Score", "G2M.Score"),
  cellcycleGenes = NULL, species = "mouse", project = "scRNAseq",
  meta = NULL, seed.use  = 10101){

  message("-------------------------")
  message("Seurat scRNAseq pipeline")
  message("-------------------------")
  message(paste("Processing:", samplename))

  # setup output dirs
  plot.dir <- file.path(outdir, "plots")
  dir.create(plot.dir, showWarnings = F)

  qc.dir <- file.path(outdir, "qc")
  dir.create(qc.dir, showWarnings = F)

  processed.counts.dir <- file.path(outdir, "processed_counts")
  dir.create(processed.counts.dir, showWarnings = F)

  Seurat.obj.dir <- file.path(outdir, "Seurat_obj")
  dir.create(Seurat.obj.dir, showWarnings = F)

  message("Importing count data")
  counts <- Read10X(counts.dir, strip.suffix = T)
  sample <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = project, meta.data = meta)

  #-----------------------
  # Sample QC
  #-----------------------
  message("Performing sample QC")
  mito.pattern <- "^mt-"
  if(species == "human" | species == "Human"){
    mito.pattern <- "^MT-"
  }

  # calculate mito percentage
  message("Calculating mito %")
  mito.genes <- grep(pattern = mito.pattern, x = rownames(x = sample@assays$RNA), value = TRUE)
  percent.mito <- Matrix::colSums(sample@assays$RNA[mito.genes, ])/Matrix::colSums(sample@assays$RNA)*100
  sample <- AddMetaData(object = sample, metadata = percent.mito, col.name = "percent.mito")

  # plot QC metrics
  message("Plotting QC metrics - Pre-filter")
  p <- VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T)
  ggsave(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_pre-filter.pdf", sep = '')), width = 9, height = 7)

  # save unfiltered obj
  # saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_01.rds", sep = "")))

  # subset based on mito count and feature cutoffs
  message("Filtering low quality cells")
  message("Starting dimensions")
  print(dim(sample))
  sample <- subset(sample, cells = colnames(sample)[which(sample$percent.mito < mitoCutoff)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA > countsLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nCount_RNA < countsUpper)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA > featuresLower)])
  sample <- subset(sample, cells = colnames(sample)[which(sample$nFeature_RNA < featuresUpper)])

  message("Filtered dimensions")
  print(dim(sample))

  # save filtered obj
  #saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_02.rds", sep = "")))

  # Feature Scatter plots
  plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mito")
  plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p <- CombinePlots(plots = list(plot1, plot2))
  ggsave(file.path(qc.dir,paste(samplename, "_count_mito_feature_QCscatter.pdf", sep = '')), width = 9, height = 7)

  p <- VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, log = T)
  ggsave(file.path(qc.dir,paste(samplename, "_10X_basic_metrics_post-filter.pdf", sep = '')), width = 9, height = 7)

  #-----------------------
  # Normalisation #1
  #-----------------------
  message(paste("Normalising counts using:", norm))
  # SCTransform workflow
  message(paste("Regressing following variables:", "percent.mito"))

  if(norm == "sctransform") {
    sample <- SCTransform(sample, vars.to.regress = "percent.mito", verbose = T, conserve.memory = F, return.only.var.genes = F)
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
  # run PCA
  sample <- RunPCA(sample, npcs = 100)
  p <- DimPlot(sample, reduction = "pca", group.by = 'Phase', label = T, label.size = 5)
  ggsave(file.path(plot.dir,paste(samplename, "_PCA_no_cc_regression.pdf", sep = '')), width = 9, height = 7)

  # Run tSNE & UMAP
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap.1.",dims,".dim", sep=''))
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne.1.",dims,".dim", sep=''))
  sample <- FindNeighbors(object = sample, verbose = T, dims = dims)

  # run a sequence of cluster resolutions
  for (i in 1:length(resolution)){
    sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
  }

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
  if(!is.null(cellcycleGenes)){
    sample <- CellCycleScoring(sample, s.features = as.vector(cellcycleGenes[[1]]), g2m.features = cellcycleGenes[[2]], set.ident = TRUE)

    # add in cell cycle difference scoring
    sample$CC.Difference <- sample$S.Score - sample$G2M.Score

    # Visualize the distribution of cell cycle markers across clusters
    p <- RidgePlot(sample, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
    ggsave(file.path(plot.dir,paste(samplename, "_ridgeplot_cell-cycle-genes.pdf", sep = '')), width = 9, height = 7)

    # run PCA on cell cycle components
    sample <- RunPCA(sample, features = c(s.genes, g2m.genes))
    p <- DimPlot(sample, reduction = "pca", group.by = 'Phase', label = T, label.size = 5)
    ggsave(file.path(plot.dir,paste(samplename, "_PCA_cell_cycle_anno.pdf", sep = '')), width = 9, height = 7)

  }

  #-----------------------
  # Normalisation #2
  #-----------------------
  message(paste("Normalising counts using:", norm))
  # SCTransform workflow
  message(paste("Regressing following variables:", regressVars))

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

  ElbowPlot(sample, ndims = 50) +
    geom_vline(aes(xintercept = elbow.dim), color = "red") +
    labs(title = paste("Elbow plot -",samplename), subtitle = paste("Predicted elbow PC =", elbow.dim))
  ggsave(file.path(plot.dir,paste(samplename,"_elbow_plot.pdf", sep = '')), width = 8, height = 6)

  # Run tSNE & UMAP on 1 to elbow dimensions
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name=paste("umap.1.",elbow.dim,".elbow", sep=''))
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = 1:elbow.dim, reduction.name=paste("tsne.1.",elbow.dim,".elbow", sep=''))

  # find neighbors based on 1 to elbow dimensions
  sample <- FindNeighbors(object = sample, verbose = T, dims = 1:elbow.dim)

  # Run tSNE & UMAP on 1 to dims
  sample <- RunUMAP(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("umap.1.",dims,".dim", sep=''))
  sample <- RunTSNE(object = sample, verbose = T, seed.use = seed.use, dims = dims, reduction.name=paste("tsne.1.",dims,".dim", sep=''))

  # run a sequence of cluster resolutions for elbow dim
  for (i in 1:length(resolution)){
    sample <- FindClusters(object = sample, verbose = T, resolution = resolution[i], random.seed = seed.use)
  }

  # save normalised and clustered dataset
  saveRDS(sample, file = file.path(Seurat.obj.dir, paste(samplename, "_cc_regressed_normalised.rds", sep = "")))

  # PCA plot reduction on clusters
  p <- DimPlot(object = sample, reduction = "pca", pt.size = 0.1, label = T, label.size = 5, dims = c(head(dims,1),tail(dims,1))) +
    NoLegend() +
    ggtitle(paste('PCA -', samplename))
  p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
  ggsave(file.path(plot.dir,paste(samplename, "_PCA.pdf", sep = '')), width = 9, height = 7)

  # UMAP reduction on clusters
  p <- DimPlot(object = sample, reduction = "umap", pt.size = 0.1, label = T, label.size = 5) +
    NoLegend() +
    ggtitle(paste('UMAP -', samplename))
  p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
  ggsave(file.path(plot.dir,paste(samplename, "_UMAP.pdf", sep = '')), width = 9, height = 7)

  # t-SNE reduction on clusters
  p <- DimPlot(object = sample, reduction = "tsne", pt.size = 0.1, label = T, label.size = 5) +
    NoLegend() +
    ggtitle(paste('t-SNE -', samplename))
  p <- AugmentPlot(p, dpi = 300, width = 8, height = 6)
  ggsave(file.path(plot.dir,paste(samplename, "_t-SNE.pdf", sep = '')), width = 9, height = 7)

  message(paste("Finished", samplename, "analysis. Enjoy your data!"))
  return(sample)
}
