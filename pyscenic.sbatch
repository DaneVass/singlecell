#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="SCENIC"
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mem=64G
#SBATCH --output="./logs/scenic-%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/scenic-%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod

pyscenic="/home/dvassiliadis/my_scripts/run_pyscenic.py"

python $pyscenic -t 24 -m 16 -i $1 -o $2 -g $3 # 24 cpus 8GB mem per cpu





