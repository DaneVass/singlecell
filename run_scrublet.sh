#!/bin/bash

logdir='./logs'

if [ ! -d $logdir ]; then
    mkdir $logdir
fi

SAMPLES=$1
doubletrate=$2

conda activate scrublet

while read line; do
    sample=$(echo $line | cut -d '/' -f9)
    
    # setup output dir
    outdir="../scrublet/${sample}-scrublet"
    if [ ! -d $outdir ]; then
        mkdir $outdir
    fi
    # toggle this to enable slurm launcher
    #cmd="sbatch scrublet.sbatch $line $outdir $sample"
    cmd="python run_scrublet.py --indir ${line} -o ${outdir} -n ${sample} -d $doubletrate
    echo $cmd
    #$cmd
done < $SAMPLES
