#!/bin/bash

#SBATCH --ntasks=40
#SBATCH --job-name="cr-aggr"
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mem=150G
#SBATCH --output="./logs/cellranger-aggr-%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/cellranger-aggr-%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod

#fastqdir=$1
#samp=$2
#outdir=$3

#transcriptome='/home/dvassiliadis/software/cellranger-3.0.1/reference/refdata-cellranger-mm10-3.0.0'
localcores=40
localmem=150

CSV="/dawson_genomics/Projects/Barcode/KF09_Clonal_competition/scripts/MLL_T0_count.csv"

cellranger aggr --id=MLL_T0 \
    --csv=$CSV \
    --nosecondary \
    --disable-ui \
    --localcores=${localcores} \
    --localmem=${localmem} \
    --normalize=mapped

$cmd
