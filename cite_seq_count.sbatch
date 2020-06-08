#!/bin/bash

#### SLURM Parameters: ####
#SBATCH --ntasks=8
#SBATCH --job-name=cite-seq-count
#SBATCH --time=00-8:00:00
#SBATCH --mem=64GB
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --output='./logs/CITE-seq_count.%j.out'
#SBATCH --error='./logs/CITE-seq_count.%j.err'
#SBATCH --partition=prod

#### Load default modules ####
module purge

#### RUN Parameters: ####

if [ ! -d ./logs ]; then
    mkdir ./logs
fi

READ1=$1
READ2=$2
SAMPLE=$3
#WL=$4

TAGS="../raw_data/HTO/antibody_tags.csv"
OUT="../raw_data/HTO/${SAMPLE}"
UNKNOWNS="unknowns.txt"
CELLS=22000
#WL="/home/dvassiliadis/projects/katie.fennell/KF05/raw_data/10X_expression/KD-COMBO-32/PRJ190741_LPRJ190741/outs/filtered_feature_bc_matrix/barcodes.tsv"

CITE-seq-Count -T 8 -R1 $READ1 -R2 $READ2 -t $TAGS -o $OUT -cells $CELLS \
    -u $UNKNOWNS --dense \
    -cbf 1 -cbl 16 -umif 17 -umil 26 \
    --max-errors 2 
#-wl $WL
