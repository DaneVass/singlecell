#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=seurat_pipeline
#SBATCH --time=8:00:00
#SBATCH --mem=128GB
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --output='seurat_pipeline-%j.out'
#SBATCH --error='seurat_pipeline-%j.err'
#SBATCH --partition='prod'

#### Load default modules ####
module purge
module load pmc-utils
module load pmc-scripts
module load R

#### Run parameters: ####
Rscript seurat_analysis_shell.R -h

Rscript seurat_analysis_shell.R -i ~/Dropbox/PostDoc_PeterMac/Projects/Barcoding_project/pbmc_1k/filtered_data \
    -o ./test \
    -n 'pbmc-test' \
    -s 'human' \
    --project 'pipeline-test' \
    -m 10 \
    -d 50 \
    --regressvars "S.Score,G2M.Score,percent.mito" \
    --seed 10101
