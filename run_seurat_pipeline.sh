#!/bin/bash

if [ ! -d ./logs ]; then
    mkdir ./logs
fi

while read line; do
    indir=$(echo $line | cut -d ',' -f 1)
    outdir=$(echo $line | cut -d ',' -f 2)
    samplename=$(echo $line | cut -d ',' -f 3)
    project=$(echo $line | cut -d ',' -f 4)
    
    # change seurat pipeline options in the sbatch
    cmd="sbatch seurat_pipeline.sbatch $indir $outdir $samplename $project"
    echo $cmd
    #$cmd
done < $1 #input is csv file containing cellranger_output path,outdir,samplename,project
