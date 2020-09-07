#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name="cellSNP"
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dane.vassiliadis@petermac.org
#SBATCH --mem=16G
#SBATCH --output="./logs/genotype_ptprc-%j.out"  # lets call the file "<jobID>"
#SBATCH --error="./logs/genotype_ptprc-%j.err"  # lets call the file "<jobID>"
#SBATCH --partition=prod_short

# Dane Vassiliadis
# PeterMac Bioinformatics 
# 20-09-07

# This script takes a BAM file from a scRNAseq 10X genomics run
# and uses cellSNP to extract and pileup SNPs of interest
# originally written to genotype cells at Ptprca to determine 
# 

source ~/software/anaconda3/etc/profile.d/conda.sh
conda activate cellSNP

BAM=$1
SAMPLENAME=$2
SNPS=$3
CELLS=$4

# run cellSNP to extract a set of known SNPs from Ptprc
#cellSNP -h
cellSNP -s $BAM -o ${SAMPLENAME}_cellSNP.vcf -R $SNPS -b $CELLS

# import to filter VCF output and generate a list of cell ID and genotype
# for import into Seurat objects as sample metadata
conda activate r-environment

#~/my_scripts/singlecell/filter_cellSNP.R -h
~/my_scripts/singlecell/filter_cellSNP.R --infile ${SAMPLENAME}_cellSNP.vcf.gz --outdir . -n ${SAMPLENAME}

echo Ptprc Genotyping complete

