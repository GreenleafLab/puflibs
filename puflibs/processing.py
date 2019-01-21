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

def combine_filenames_split(filenames, avoid_elements=['ann', 'filt', 'bedGraph', 'tracks', 'txt', 'gz', 'counts', 'pkl']):
    """Combine a set of filenames, removing any duplicate elements separated by periods"""

    name_list = list(itertools.chain(*[os.path.basename(filename).split('.') for filename in filenames]))
    _, idx = np.unique(name_list, return_index=True)
    name_unique = [name_list[i] for i in np.sort(idx)]
    basename = '.'.join([s for s in name_unique if s not in avoid_elements])

    return basename

def load_motif(filename):
    """Load a homer motif"""
    with open(filename) as f:
        lines = f.read_lines()
    
    header = lines[0][1:].split()
    

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

def find_roc_data(data, prediction_column, actual_column, predicted_up=True,):
    """data is a dataframe containing the prediction column (should be quantitative variable) and actual_column (should be bool)"""

    data_sub = data.dropna(subset=[prediction_column, actual_column]).copy()
    threshold_values = data_sub.loc[:, prediction_column].quantile(np.linspace(0, 1, 500)).dropna().unique()
    if len(threshold_values)<=2:
        print "data.loc[:, '%s'] had only two or fewer unique values"%prediction_column
        return

    actual_positive = data_sub.loc[:, actual_column]

    roc_curve = {}
    for val in threshold_values:
        #num_false_positives, num_true_positives = pd.concat([expression_data.sig_up, expression_data.occupancy >= val], axis=1).loc[~expression_data.occupancy.isnull()].groupby([ 'occupancy', 'sig_up']).size().loc[True]
        if predicted_up:
            predicted_positive = data_sub.loc[:, prediction_column] >= val
        else:
            predicted_positive = data_sub.loc[:, prediction_column] <= val
        tpr, fpr = get_tpr_fpr(predicted_positive, actual_positive)

        roc_curve[val] = pd.Series({'tpr':tpr, 'fpr':fpr})
    roc_curve = pd.concat(roc_curve).unstack()
    return roc_curve

def get_tpr_fpr(predicted_bool_vec, actual_bool_vec):
    """Given two vectors, one with True/False predicted values, one with True/Fals,e actual values, and return tpr and fpr"""
    actual_bool_vec = actual_bool_vec.astype(bool)
    predicted_bool_vec = predicted_bool_vec.astype(bool)
    num_actual_positive = actual_bool_vec.sum()
    num_actual_negative = actual_bool_vec.shape[0] - num_actual_positive
    num_false_positives = (predicted_bool_vec&(~actual_bool_vec)).sum()
    num_true_positives = (predicted_bool_vec&actual_bool_vec).sum()    
    return float(num_true_positives)/num_actual_positive, float(num_false_positives)/num_actual_negative


def get_counts_from_counts_table(data_table, interval_radius=40, offset=15):
    """For a data table, with each row giving the count per base pair of that site,
    count the total between the start and end base pair spanning 2*interval_radius, offset by offset"""
    window_size = data_table.shape[1]
    start_loc = int(window_size/2-interval_radius-offset)
    end_loc = start_loc + interval_radius*2
    return data_table.iloc[:, start_loc:end_loc].sum(axis=1)

        
def trapz(y, x):
    """Return the trapezoidal integrated input"""
    total = 0
    for i in range(1, len(y)):
        total += (y[i] + y[i-1])/2.*(x[i] - x[i-1])
    return total