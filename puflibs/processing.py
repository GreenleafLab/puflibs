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
from puflibs import variables

def load_bed(filename, additional_cols=None):
    """load bed. ASSUMES NO HEADER"""
    cols = variables.bed_fields
    if additional_cols:
        cols = cols + additional_cols

    regions = pd.read_table(filename, header=None, names=cols)
    col_names = regions.columns.tolist()
    if col_names[0][0] == '#':
        col_names = [col_names[0][1:]] + col_names[1:]
    regions.columns = col_names
    return regions

def load_bed_header(filename):
    """Need to remove a prepended #. ASSUMES HEADER"""
    regions = pd.read_table(filename)
    col_names = regions.columns.tolist()
    if col_names[0][0] == '#':
        col_names = [col_names[0][1:]] + col_names[1:]
    regions.columns = col_names
    return regions

def save_bed(regions, filename):
    """Need to add a prepended #"""
    regions.to_csv(filename, sep='\t', header=False, index=False)
    
def save_bed_header(regions, filename):
    """Need to add a prepended #"""
    regions = regions.copy()
    col_names = regions.columns.tolist()
    if col_names[0][0] != '#':
        col_names = ['#' + col_names[0]] + col_names[1:]
    regions.columns = col_names
    regions.to_csv(filename, sep='\t', index=False)

def combine_filenames(filenames):
    """Combine a set of filenames, removing any prepending str from later filenames common with first filename"""
    basenames = [os.path.basename(os.path.splitext(filename)[0]) for filename in filenames]
    ref_basename = basenames[0]
    basenames_filtered = [ref_basename]
    for i in range(1, len(basenames)):
        # find the parts of the current basename that match the previous
        curr_basename = basenames[i]
        common_str = ''
        for c, c_ref in zip(curr_basename, ref_basename):
            if c!=c_ref:
                break
            common_str += c
            
        new_basename = curr_basename[len(common_str):]
        basenames_filtered.append(new_basename)
    
    return basenames_filtered

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


    
