#!/home/dvassiliadis/software/anaconda3/envs/pyscenic/bin/python
# run_pyscenic

# Dane Vassiliadis
# Peter MacCallum Cancer Center
# June 2020

# automates scenic run for python implementation
# based on the vignette on the py_scenic github

#------ Initialise env
import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from dask.distributed import Client, LocalCluster

import seaborn as sns
from optparse import OptionParser
import sys

def main():
    # parse input args
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Input counts matrix file")
    parser.add_option("-o", "--output", dest="output", help="Output file")
    parser.add_option("-d", "--tmpdir", dest="tmpdir", help="Temporary working dir", default=".")
    parser.add_option("-t", "--threads",dest="threads", help="Number of processes to spawn", type="int", default=8)
    parser.add_option("-g", "--genome", dest="genome", help="Reference genome (mm10 or hg19)", default="hg19")
    parser.add_option("-m", "--memory", dest="memory", help="Memory limit (in GB)", default=8, type="float")
    (options, args) = parser.parse_args()

    if options.input is None or options.output is None:
        sys.exit("Input and output files must be specified.")

    DATA_FOLDER=options.tmpdir
    RESOURCES_FOLDER="/home/dvassiliadis/software/pySCENIC/resources"
    DATABASE_FOLDER = "/home/dvassiliadis/software/pySCENIC/resources"
    DATABASES_GLOB = os.path.join(DATABASE_FOLDER, options.genome + "*.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, ("motifs-v9-nr.hgnc-m0.001-o0.0.tbl" if options.genome == "hg19" else "motifs-v9-nr.mgi    -m0.001-o0.0.tbl"))
    TFS_FNAME = os.path.join(RESOURCES_FOLDER, ('hs_hgnc_tfs.txt' if options.genome == "hg19" else "mm_mgi_tfs.txt"))
    SC_EXP_FNAME = options.input
    REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
    MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")

    # import the expression matrix and transpose it to have genes as columns and cells as rows
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
    print(ex_matrix)
    ex_matrix.shape

    tf_names = load_tf_names(TFS_FNAME)
    print(tf_names)
    db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    dbs

    print("Running Phase I - Calculate cellular adjacencies")
    cluster = LocalCluster(n_workers = options.threads, threads_per_worker=1, memory_limit=options.memory * 10e9)
    client = Client(cluster)
    
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True, client_or_address = client)
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    
    # Phase II: Prune modules for targets with cis regulatory footprints (aka RcisTarget)
    print("Running Phase II - Prune modules for targets with cis regulatory footprints (aka RcisTarget)")
    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address = client)
    
    # Create regulons from this table of enriched motifs.
    print("Generating regulons from enriched motifs")
    regulons = df2regulons(df)
    
    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

    # The clusters can be leveraged via the dask framework:   
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address = client)

    df = load_motifs(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "rb") as f:
        regulons = pickle.load(f)

    # Phase III: Cellular regulon enrichment matrix (aka AUCell)
    print("Running Phase III - Regulon enrichment using AUCell")
    auc_mtx = aucell(ex_matrix, regulons, num_workers=options.threads)
    #sns.clustermap(auc_mtx, figsize=(8,8))
    auc_mtx.to_csv(options.output)

if __name__ == "__main__":
    main()
