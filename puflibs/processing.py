import os
import sys
import time
import re
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import scipy.stats as st
import warnings
import collections
import pickle
import seaborn as sns
import itertools


def process_go_term_gene_string(s):
    """Process the go term gene string."""
    snew = s.strip(']"').lstrip('"[')[:]
    genes = [a[1:a.find('  -  ')] for a in snew.split(',') if a.find('  -  ')>=0]
    return genes

def load_go_file_with_errors(filename):
    """This was designed to load the outpt of gorilla go terms. especially if there are things with weird line breaks."""
    with open(filename) as f:
        line = f.readlines()
    if len(line)>1:
        print '%d lines in file.'%len(line)
        line = '\r'.join(line)
    elif len(line)==1:
        line = line[0]
    lines = line.split('\r')
    header = ['_'.join(s.lower().replace('-', '').split()) for s in lines[0].split('\t')]
    mat = {}
    genes = {}
    counter = 0
    for line in lines[1:]:
        vec = line.split('\t')
        if vec[0].split(':')[0] == 'GO':
            mat[counter] = pd.Series(vec[:-1], index=header[:-1])
            genes[counter] = process_go_term_gene_string(vec[-1])
            counter += 1
        else:
            genes[counter-1] = genes[counter-1] + process_go_term_gene_string(vec[0])

    mat = pd.concat(mat).unstack()
    mat.loc[:, 'genes'] = pd.Series(genes)
    
    return  mat.set_index('go_term')

def normalize_clip_signal_tpm(data, num_reads_total, count_cols=['input', 'rep1', 'rep2']):
    """Using the count data and tpm data, return normalize fold enrichment per site."""
    num_reads_total = pd.read_pickle('analysis/output/total_bamcounts.pkl')
    weights_sample = 1./num_reads_total*num_reads_total.mean()
    weights_site = 1/data.tpm*data.tpm.mean()
    
    data_norm = {}
    for col in count_cols:
        data_norm[col + '_norm'] = data.loc[:, col]*weights_site*weights_sample.loc[col]
    return pd.concat(data_norm, axis=1,)
    
