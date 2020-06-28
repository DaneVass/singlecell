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
module load R/3.6.0

#### Run parameters: ####
Rscript seurat_analysis_shell.R -h

Rscript seurat_analysis_shell.R -i $1 \
    -o $2 \
    -n $3 \
    -s 'mouse' \
    --project $4 \
    -m 10 \
    -d 50 \
    --regressvars "S.Score,G2M.Score,percent.mito" \
    --seed 10101
