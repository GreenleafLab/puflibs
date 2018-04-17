#!/usr/bin/env python
""" Make figures for paper.

Sarah Denny """

##### IMPORT #####
import numpy as np
import pandas as pd
import os
import argparse
import pickle
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess

### MAIN ###
#set up command line argument parser
parser = argparse.ArgumentParser(description='bootstrap on off rate fits')
parser.add_argument('--data', default="paper/01_expt/filename_table.dat",
                   help='file giving filenames')
parser.add_argument('--key', nargs='*',
                   help='key of experiment to look at')
parser.add_argument('--mode', 
                   help='what to do with it')
parser.add_argument('--out', 
                   help='indication of output file')
parser.add_argument('-f', '--force_overwrite', action="store_true", 
                   help='whether to overwrite output')

if __name__ == '__main__':
    args = parser.parse_args()
    
    filenames = ['analysis/output/split_%d/all_unprocessed_st_merged.0%d.hPUM2_all.random_1e+06.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.combined_data.gz'%(i, i) for i in range(10)]
    
    all_data = pd.concat({i:pd.read_table(filename, compression='gzip', index_col=0).reset_index().rename(columns={'index':'name'}) for i, filename in enumerate(filenames)})
    all_data.loc[:, 'name'] = all_data.name + '_' + pd.Series({idx:'%d'%idx[0] for idx in all_data.index.tolist()})
    all_data_comp = all_data.sort_values('score').groupby(['chrm', 'start', 'stop']).first()
    
    all_data_comp.reset_index().to_csv('analysis/output/all_unprocessed_st_merged.hPUM2_all.random_1e+06.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.combined_data.gz', compression='gzip', index=False)
    
    
    # group by
    filename = 'analysis/clip/split_7/all_unprocessed_st_merged.07.hPUM2_all.random_1e+06.07.ann.filt.input.ENCFF786ZZB.R2.bedGraph.500.tracks.txt.gz'
    vec = []
    counts = pd.concat([chunk.iloc[:, start_loc:end_loc].sum(axis=1) for chunk in pd.read_csv(filename, compression='gzip', index_col=0, chunksize=10000)])

               
    
    