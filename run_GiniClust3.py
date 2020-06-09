""" 
GiniClust3.py

Dane Vassiliadis. Peter Mac. Mar 2020

Run GiniCLust3 on a single cell counts matrix

"""

#------------------
# Load packages
#------------------

import re
import sys
import time
import os
import scanpy as sc
import numpy as np
import giniclust3 as gc
import anndata
import gzip
import shutil
from optparse import OptionParser

startTime = datetime.now()


#------------------
# Define options
#------------------
usage = "USAGE: GiniCLust3.py -i [input-matrix] -o [output-dir]"

parser = OptionParser()
parser.add_option('-i', '--input', type = 'string', dest = 'input', help="Path to input fastq file.")
parser.add_option('-o', '--output', type = 'string', dest = 'outfile', help="Path to output fastq file.")
#parser.add_option('--upconstant', type = 'string', dest = 'upconstant', help="Upstream (i.e. 5') constant region.")

(options,args) = parser.parse_args()

#------------------
# Parse inputs
#------------------


def 

# Read single cell file
adataRaw=sc.read_csv("./data/GSM1599495_ES_d0_biorep_techrep1.csv",first_column_names=True)
Filter expression matrix

sc.pp.filter_cells(adataRaw,min_genes=3)
sc.pp.filter_genes(adataRaw,min_cells=200)
Format expression matrix

# example csv file is col:cells X row:genes. Skip this step if the input matrix is col:genes X row:cells
adataSC=anndata.AnnData(X=adataRaw.X.T,obs=adataRaw.var,var=adataRaw.obs)

# Normalization
sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)
Perform GiniIndexClust

gc.gini.calGini(adataSC) ###Calculate Gini Index
adataGini=gc.gini.clusterGini(adataSC,neighbors=3) ###Cluster based on Gini Index

# Perform FanoFactorClust
gc.fano.calFano(adataSC) ###Calculate Fano factor
adataFano=gc.fano.clusterFano(adataSC) ###Cluster based on Fano factor

# ConsensusClust
consensusCluster={}
consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
gc.consensus.generateMtilde(consensusCluster) ###Generate consensus matrix
gc.consensus.clusterMtilde(consensusCluster) ###Cluster consensus matrix
np.savetxt("final.txt",consensusCluster['finalCluster'], delimiter="\t",fmt='%s')

# UMAP visualization
adataGini.obs['final']=consensusCluster['finalCluster']
adataFano.obs['final']=consensusCluster['finalCluster']
gc.plot.plotGini(adataGini)
gc.plot.plotFano(adataFano)