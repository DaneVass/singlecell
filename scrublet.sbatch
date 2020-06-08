#!/bin/bash

#### SLURM Parameters: ####
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --job-name=scrublet
#SBATCH --time=00-2:00:00
#SBATCH --mem=32GB
#SBATCH --mail-user='dane.vassiliadis@petermac.org'
#SBATCH --mail-type=ALL
#SBATCH --output='logs/scrublet-%j.out'
#BATCH --error='logs/scrublet-%j.err'
#SBATCH --partition=prod_short

#### Load default modules ####
module purge
#conda activate scrublet

#### RUN Parameters: ####

INPUT=$1
OUTDIR=$2
SAMPLE=$3

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

python run_scrublet.py -i $INPUT -o $OUTDIR -n $SAMPLE -d '0.1'

