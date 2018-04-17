import os
import subprocess
import logging
import itertools
import pandas as pd
import numpy as np
import ipdb
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from puflibs import variables
from fittinglibs import seqfun
import cPickle as pickle

from joblib import Parallel, delayed





def process_bed(bed_data, ncores):
    """For the data in bed data, group by the name."""
    bed_data = bed_data.copy()
    index = pd.MultiIndex.from_tuples([s.split('/') for s in bed_data.name])
    bed_data.index = index
    names_all = bed_data.index.levels[0].tolist()
    names_split = [vec.tolist() for vec in np.array_split(names_all, ncores)]
    results = Parallel(n_jobs=ncores, verbose=10)(delayed(process_set)(bed_data, names) for names in names_split)
    
    processed = pd.concat([df[0] for df in results])
    not_processed = pd.concat([df[1] for df in results])
     
    return processed, not_processed


def process_set(bed_data, names):
    """For a set of names, return the processed and unprocessed results."""
    not_processed = []
    processed = []
    for name in names:
        group = bed_data.loc[name].copy()
        # if there was only one read, or if there is more than one chromosome, don't process
        if len(group)==1 or len(group.chrm.unique())!=1 or len(group)>2:
            not_processed.append(group)
        else:
            chrm, _, _, _, score, strand = group.loc['1']
            if strand == '+':
                start = group.loc['1'].start
                stop = group.loc['2'].stop
            elif strand == '-':
                start = group.loc['2'].start
                stop = group.loc['1'].stop
            processed.append(pd.Series([chrm, start, stop, name, score, strand], index=bed_data.columns))
    if not processed:
        processed = pd.Series(index=bed_data.columns)
    else:
        processed = pd.concat(processed, axis=1).transpose()
    if not not_processed:
        not_processed = pd.Series(index=bed_data.columns)
    else:
        not_processed = pd.concat(not_processed).reset_index(drop=True)
    return processed, not_processed

def test_process(bed_data):
    """this is a function"""
    
    return {}

if __name__=="__main__":
    # import args
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input bed file, with name column giving the read')
    parser.add_argument('-o', '--output', help='output bed file, merged for the two reads')
    parser.add_argument('-n', '--ncores', help='number of cores', type=int, default=10)



    args = parser.parse_args()
    
    chunksize = 1E6
    bedfields = ['chrm', 'start', 'stop', 'name', 'score', 'strand']
    
    ncores = args.ncores
    """
    # make sure files don't already exist, or you'll append on to previous older version of the file
    filename1 = args.output + '.processed.bed'
    filename2 = args.output + '.unprocessed.bed'
    for filename in [filename1, filename2]:
        if os.path.exists(filename):
            os.remove(filename)
    
    # intiate output files
    f1 = open(args.output + '.processed.bed', 'a')
    f2 = open(args.output + '.unprocessed.bed', 'a')
    
    # go through chunks, group and write to output
    #for chunk in pd.read_table(args.input, names=bedfields, header=None, chunksize=chunksize):
    """
    chunk =  pd.read_table(args.input, names=bedfields, header=None)
    #pickle.dump(test_process, 'test.p')
    results = Parallel(n_jobs=ncores,)(delayed(test_process)(row) for idx, row in chunk.iterrows())
    ipdb.set_trace()
    """
        bed_data = chunk.copy()
        bed_data.index = pd.MultiIndex.from_tuples([s.split('/') for s in bed_data.name])
        names_all = bed_data.index.levels[0].tolist()
        names_split = [vec.tolist() for vec in np.array_split(names_all, ncores)]
        results = Parallel(n_jobs=ncores, verbose=10)(delayed(test_process)(bed_data, names) for names in names_split)
        processed = pd.concat([df[0] for df in results])
        not_processed = pd.concat([df[1] for df in results])
        
       
        

        processed.to_csv(f1, header=False, index=False, sep='\t')
        not_processed.to_csv(f2, header=False, index=False, sep='\t')

    f1.close()
    f2.close()
    """    