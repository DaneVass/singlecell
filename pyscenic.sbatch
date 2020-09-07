#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name="SCENIC"
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mem=128G
#SBATCH --output="./logs/scenic-%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/scenic-%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod

pyscenic="/home/dvassiliadis/my_scripts/run_pyscenic.py"

python $pyscenic -t 12 -m 10 -i $1 -o $2 -g $3 -d $4 # 24 cpus 8GB mem per cpu





