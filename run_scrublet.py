#!~/software/anaconda3/envs/scrublet/bin/python

# Dane Vassiladis 
# 19-10-24

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import gzip
from optparse import OptionParser
import subprocess

usage = """ run scrublet script

USAGE: python run_scrublet.py [OPTIONS]

Arguments:
-i indir:       Path to input directory. Should contain raw counts.mtx and barcodes.csv
-o outdir:      Path to output directory. Results will be written here (default = './')
-n samplename:  Desired prefix to name output files.
-d doubletrate: Expected rate of doublets given technology and cells loaded. Default = 0.06 (6%)
"""

# setup argument options
parser = OptionParser()
parser.add_option('-i', '--indir', type = 'string', dest = 'indir', help='Path to input directory. Should contain raw counts.mtx and barcodes.csv')
parser.add_option('-o','--outdir', type = 'string', dest = 'outdir', default = "./", help = "Path to output directory. Results will be written here (default = './')")
parser.add_option('-n', '--samplename', type = 'string', dest = 'samplename', help = "sample name for labelling output files", default = 'scrublet')
parser.add_option('-d', '--doubletrate', type = 'string', dest = 'doubletrate', help = "expected doublet rate based on the platform used and the number of cells loaded into the experiment.")
(options,args)=parser.parse_args()

# check inputs
if options.indir is not None:
    indir = options.indir
    try:
        os.path.isdir(indir)
    except:
        print("Please set path to an existing directory containing matrix, barcode, and UMI counts files. Exiting")
        sys.exit()

# check output dir
if options.outdir is not None:
    outdir = options.outdir
    try:
        os.path.isdir(outdir)
    except:
        print("An output directory does not exist at {}. Creating".format(outdir))
        os.mkdir(outdir)
else:
    outdir = "./"
    print("No output directory given. Defaulting to {}".format(os.getcwd()))

# Set doublet rate
if options.doubletrate is not None:
    doubletrate = float(options.doubletrate)
    print("Expected doublet rate has been set to {}".format(doubletrate))
else:
    doubletrate = 0.06
    print("No doublet rate given, defaulting to 0.06 (6%)")

# Setup prefix
prefix = options.samplename

# Matplotlib options
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

# import counts data
print("Analysing {} with Scrublet".format(prefix))
counts_matrix = scipy.io.mmread(indir + '/matrix.mtx.gz').T.tocsc()

# decompress the features.tsv file first
features = indir + '/features.tsv.gz'
subprocess.call(['gzip', '-d', features])
genes = np.array(scr.load_genes(indir + '/features.tsv', delimiter='\t', column=2))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

# Setup Scrublet object
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate = doubletrate)

# Generate doublet scores
doublet_scores, predicted_doublets = scrub.scrub_doublets()
#doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
#                                                          min_cells=3, 
#                                                          min_gene_variability_pctl=85, 
#                                                          n_prin_comps=30)

# plot histogram of simulated vs observed doublets
scrub.plot_histogram()
plt.savefig(os.path.join(outdir, prefix + "_scrublet_histogram.pdf"))

# Run dimensinal reductions
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

# # Uncomment to run tSNE - slow
#print('Running tSNE...')
#scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))

# # Uncomment to run force layout - slow
#print('Running ForceAtlas2...')
#scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5., n_iter=1000))

# Plot dimensional reductions
scrub.plot_embedding('UMAP', order_points=True);
plt.savefig(os.path.join(outdir, prefix + "_scrublet_UMAP.pdf"))

#scrub.plot_embedding('tSNE', order_points=True);
#plt.savefig(os.path.join(outdir, prefix + "_scrublet_tSNE.pdf"))

#scrub.plot_embedding('FA', order_points=True);
#plt.savefig(os.path.join(outdir, prefix + "_scrublet_ForceLayout.pdf"))

# Write predicted doublets out to a file
np.savetxt(os.path.join(outdir, prefix + "_scrublet_predicted_doublets.csv"), predicted_doublets, delimiter=",")
np.savetxt(os.path.join(outdir, prefix + "_scrublet_doublet-scores.csv"), doublet_scores, delimiter=",")

print("Scrublet analysis of {} complete!".format(prefix))
