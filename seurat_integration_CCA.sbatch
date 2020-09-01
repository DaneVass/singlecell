#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name="seu_int_CCA"
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mem=128G
#SBATCH --output="./logs/seruat-integration-cca-%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/seurat-integration-cca-%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod

seuratobj=
outdir="./"

Rscript seurat_integration.R $seuratobj $seuratnames $outdir


