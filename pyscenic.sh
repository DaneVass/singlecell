#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name="SCENIC"
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mem=200G
#SBATCH --output="./logs/scenic-%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/scenic-%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod

source ~/software/anaconda3/etc/profile.d/conda.sh
conda activate scenic

EXP_MAT="/home/dvassiliadis/my_scripts/singlecell/GSE60361_C1-3005-Expression.csv"
TF="/home/dvassiliadis/software/pySCENIC/resources/mm_mgi_tfs.txt"

# start scenic run

# Phase 1 - Run GrnBoost2
# Outputs a list of adjacencies between a TF and its targets
arboreto_with_multiprocessing.py -o adjacencies.csv \
    --seed 10101 -t --num_workers 20 -m grnboost2 ${EXP_MAT} ${TF}

# Phase 2 - Define Regulons

FEATHER_DB=$(find /home/dvassiliadis/software/pySCENIC/resources -name "mm10*.feather")
MOTIFS="/home/dvassiliadis/software/pySCENIC/resources/motifs-v9-nr.mgi-m0.001-o0.0.tbl"

pyscenic ctx adjacencies.csv ${FEATHER_DB} \
    -t \
    --annotations_fname ${MOTIFS} \
    --expression_mtx_fname ${EXP_MAT} \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 20
# in earlier versions the default behavior was to mask dropouts, 
# while in v0.9.18 and later, the correlation is performed using the entire set of cells (including those with zero expression)
# to use all cells remove --mask-dropouts

# Phase 3 - AU Cell analysis
OUTFILE="./outfile.csv"
pyscenic aucell $EXP_MAT reg.csv --output $OUTFILE --num_workers 20 --seed 10101 -t

echo SCENIC COMPLETE
