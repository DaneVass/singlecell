#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --job-name="cellranger"
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mem=100G
#SBATCH --output="./logs/%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod

fastqs=$1
outdir=$2

transcriptome="/home/mzethoven/cellranger/refdata-cellranger-hg19-1.2.0"
localcores=30
localmem=100

cellranger count --id=$outdir --fastqs=$fastqs --transcriptome=$transcriptome --localcores=$localcores --localmem=$localmem

