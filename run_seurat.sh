#!/bin/bash

Rscript seurat_analysis_shell.R -i ~/Dropbox/PostDoc_PeterMac/Projects/Barcoding_project/pbmc_1k/filtered_data \
    -o ./test \
    -n 'pbmc-test' \
    -s 'human' \
    --project 'pipeline-test' \
    -m 10 \
    -d 50 \
    --regressvars "S.Score,G2M.Score,percent.mito" \
    --seed 10101


