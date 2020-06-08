
""" 
Dane Vassiliadis. Peter Mac. Jan 2019

Map SPLINTR barcodes to single cell RNAseq cell barcodes

Takes a fastq containing reads from DNA barcoded cells, 
      a reference whitelist of SPLINTR barcodes and
      a list of cell barcodes from a 10X experiment processed with cellranger count output

Returns a file containing the cell barcode and the most likely DNA barcode.  

Assumes following barcode structure:
5' - tgaccatgtacgattgactaNNSWSNNWSWNNSWSNNWSWNNSWSNNWSWNNSWSNNWSWNNSWSNNWSWNNSWSNNWSWtgctaatgcgtactg - 3'


From a fastq input, filter out reads with a correct 5' constant region that dont contain N's
and return a fasta file containing the remaining barcode region.

From a fastq input, filter out reads that meet the following criteria:
- Exact match to the specified 5' and 3' constant regions
- [OPTIONAL] A correct sample index
- Average quality >= 30 across barcode
- No ambiguous N residues
- Barcode length at least 30

Optionally prefix barcodes with the sample index. 
Specify --index as tab delimited file with columns ID and Sequences. e.g.
ID	            Sequences
400nM-IBET-1	    ATCACG

"""
#%%
# imports
import re
import sys
import time
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from datetime import datetime
import gzip
import shutil
from optparse import OptionParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%

startTime = datetime.now()

# paths for testing
#path = "/Users/vassiliadisdane/Dropbox/PostDoc_PeterMac/Projects/Dawson_Lab/Katie/barcoding_experiment/scripts/test.fastq"

usage = "USAGE: filter_trim_barcode_reads.py -i [input-fastq] -o [output-filename] --upconstant [5' constant-region] --downconstant [3' constant region] -r [reference-barcode-whitelist] -b [10X barcodes file]"

# Define command line options
parser = OptionParser()
parser.add_option('-i', '--input', type = 'string', dest = 'input', help='Path to input fastq file')
parser.add_option('-o', '--output', type = 'string', dest = 'outfile', help='Path to output fastq file')
parser.add_option('-u', '--upconstant', type = 'string', dest = 'upconstant', help="5' constant region")
parser.add_option('-d','--downconstant', type = 'string', dest = 'downconstant', help="3' constant region")
parser.add_option('-r','--reference', type = 'string', dest = 'reference', help='Path to file containing reference whiltelist of DNA barcodes')
parser.add_option('-b', '--barcodes', type = 'string', dest = 'cellbarcodes', help='Path to cellranger counts output barcodes.tsv for this 10X run.')


(options,args) = parser.parse_args()

# parse input
if options.input is not None:
    try:
        os.path.isfile(options.input)
        path = options.input
        print("Input file set to ", os.path.basename(path))
    except:
        print("no input fastq defined. exiting")
        print(usage)
        sys.exit(1)

# parse output
if options.outfile is not None:
    print("Writing output to ", options.outfile)
    outfile_name = options.outfile
else:
    sys.exit("no output file defined. exiting")
    print(usage)  

# parse constant regions
if options.upconstant is not None:
    upstream_constant = options.upconstant
    print("5' constant region: ", upstream_constant.upper())
else:
    upstream_constant = "GATCCTGACCATGTACGATTGACTA"
    print("no 5' constant region given, setting constant region to: ", upstream_constant.upper())

if options.downconstant is not None:
    downstream_constant = options.downconstant
    print("3' constant region: ", downstream_constant.upper())
else:
    downstream_constant = "TGCTAATGCGTACTGACTAGGCTAG"
    print("no 3' constant region given, setting constant region to: ", downstream_constant.upper())

# parse reference file
if options.reference is not None:
    try:
        os.path.isfile(options.reference)
    except:
        sys.exit("no valid reference DNA barcode whitelist specified. exiting")
    print("reference DNA barcode whitelist: ", options.reference)
    
    # create list of cell barcodes
    reference = pd.read_csv(options.reference, sep=' ', header=0)   
else:
    sys.exit("no reference DNA barcode whitelist given. exiting")
    print(usage)

# parse cellranger count cell barcodes
if options.cellbarcodes is not None:
    try:
        os.path.isfile(options.cellbarcodes)
    except:
        sys.exit("")
    print("10X cell barcodes: ", options.cellbarcodes)    
    cellbarcodes = options.cellbarcodes
    # create list of cell barcodes
    cell_barcodes = open(cellbarcodes, 'r')
    cell_barcodes = cell_barcodes.readlines()
    print(cell_barcodes)
    # create pandas df
    barcodes.df = pd.DataFrame({ '10X_Cell_Barcode' : cell_barcodes})
    barcodes.df.head()
else:
    sys.exit("no list of 10X cell barcodes defined. exiting")
    print(usage)

#%%

# Parsing function
def parse_fastq(file, upstream_constant, downstream_constant):
            
    # get index of sample
    if options.index_file is not None:
        index_file = open(indexes, "r")
        for line in index_file:
            if line.split('\t')[0] == str(filename):
                sample_index = str(line.split('\t')[1])
                print("Sample index: ", sample_index)
                # make biopython seq object
                index = Seq(sample_index.strip('\n'))
    # parse constant regions              
    upstream_constant = upstream_constant.upper()
    print("5' constant region: 5'", str(upstream_constant), "3'")
    
    downstream_constant = downstream_constant.upper()
    print("3' constant region: 5'", str(downstream_constant), "3'")
    print("")
    
    # parse fastq, only use index if given
    fastq_parser = SeqIO.parse(file, "fastq")
    outfile = open(outfile_name, "w")
    
    def filter_and_trim_read(fastq):
        # collect counts of total and filtered reads
        total_count = 0
        filtered_count = 0
        
        # main stuff happens here
        for record in fastq_parser:
            total_count += 1
            # take only complete reads with correct start and end constant regions
            if str(upstream_constant) in str(record.seq):
                if str(downstream_constant) in str(record.seq):
                    #print("ID ", record.id)
                    #print("Desc ", record.description)
                    #print("Seq ", record.seq)
                    read = str(record.seq)
                    start = read.index(str(upstream_constant)) + len(upstream_constant)
                    end = read.index(str(downstream_constant))
                    trimmed_read = record[start:end]
                    #print(trimmed_read.seq)
                    if options.index_file is not None:
                        desc = str(trimmed_read.description)
                        # if the read index is correct
                        if desc.strip().split(":")[9] == sample_index.strip('\n'):
                            # prefix read with sample index
                            mutable_read = MutableSeq(trimmed_read)
                            final_read = index + mutable_read
                            # replace trimmed_read with new sequence
                            trimmed_read = SeqRecord(final_read, id = trimmed_read.id, description = trimmed_read.description)
                            # check length is OK and return read
                            if len(trimmed_read) == 68:
                                yield(trimmed_read)
                                filtered_count += 1
                    
                    # else if no index check length is OK and return barcode only
                    elif len(trimmed_read) == 60:
                        read_quals = trimmed_read.letter_annotations["phred_quality"]
                        avg_qual = sum(trimmed_read.letter_annotations["phred_quality"]) / len(trimmed_read)
                        print(avg_qual)
                        if avg_qual >= 30:
                            if "N" not in str(trimmed_read.seq):
                                yield trimmed_read
                                filtered_count += 1

        print("")
        print("Total reads parsed: ", total_count)
        print("Number of reads kept: ", filtered_count)
        print("Percentage of reads kept: ", (filtered_count / total_count) * 100)

    # write filtered and trimmed reads to file
    SeqIO.write(filter_and_trim_read(fastq_parser), outfile, "fastq")
    outfile.close()

parse_fastq(path, upstream_constant, downstream_constant)
    
# gzip the outfile, can probably do this in one step?
    #with open(outfile_name, 'rb') as f_in:
    #    with gzip.open(outfile + '.gz' , 'wb') as f_out:
    #        shutil.copyfileobj(f_in, f_out)

# finish and print runtime
print("Script runtime: ", datetime.now() - startTime)
    #%%
