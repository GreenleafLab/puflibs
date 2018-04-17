#!/usr/bin/env python
#

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

# import args
parser = argparse.ArgumentParser()
parser.add_argument('--bed', help='file giving bed file')
parser.add_argument('--out', help='directory to save to')
parser.add_argument('--mode', help='[pyatac|analysis]')
args = parser.parse_args()


# filenames for clip bedgraphs
dirname_minus = 'CLIP/hPUM2/unique_minus_strand'
minus_strand_inputs = {'rep1':'%s/rep1.ENCFF452BKL.bedGraph.gz'%dirname_minus,
                        'rep2':'%s/rep2.ENCFF657MOM.bedGraph.gz'%dirname_minus,
                        'input':'%s/input.ENCFF034QCL.bedGraph.gz'%dirname_minus}


# filenames for clip bedgraphs
dirname_plus = 'CLIP/hPUM2/unique_plus_strand'
plus_strand_inputs = {'rep1':'%s/rep1.ENCFF769SCB.bedGraph.gz'%dirname_plus,
                        'rep2':'%s/rep2.ENCFF133NMP.bedGraph.gz'%dirname_plus,
                        'input':'%s/input.ENCFF880UXZ.bedGraph.gz'%dirname_plus}

# run signal
outfile_dict = {}
for input_dict, strand in zip([plus_strand_inputs, minus_strand_inputs], ['plus_strand', 'minus_strand']):
    for key, filename in input_dict.items():
        outfile_dict[(strand, key)] = '%s/clip_signal/%s/%s.%s'%(args.out, strand, key, os.path.splitext(os.path.basename(args.bed) )[0])
        command_txt = ('nohup pyatac signal --bed %s --bg %s --sizes /shr/gSizes/hg38.genomsize '
                       '--out %s --all --up 100 --down 100 --strand 6 >> %s.out &'%
                       (args.bed, filename, outfile_dict[(strand, key)], 'nohup.%s.%s'%(strand, key)))
        if args.mode == 'pyatac':
            print command_txt
            os.system(command_txt)

# find enrichment
if args.mode=='analysis':

    data = pd.concat({name:pd.read_csv(filename + '.tracks.txt.gz', compression='gzip', header=None) for name, filename in outfile_dict.items()})

    data_agg = data.sum(axis=1).unstack(level=0).fillna(0)
    data_agg2 = (data_agg.plus_strand - data_agg.minus_strand).unstack(level=0)
    min_val = data_agg2.rep1.replace([0], np.nan).min() + data_agg2.rep2.replace([0], np.nan).min()
    min_val_input = data_agg2.input.replace([0], np.nan).min()
    clip_enrichment = (data_agg2).apply(lambda x: np.log2((x.rep1+x.rep2+min_val)/(2.*(x.input))), axis=1).replace([-np.inf, np.inf], np.nan)
    
    bed_fields = ['chrm', 'start', 'stop', 'name',  'score', 'strand',];
    bed_data = pd.read_table(args.bed, header=None, names=bed_fields + ['annotation', 'gene', 'exon'], index_col='name')
    
    clip_enrichment.index = bed_data.index
    