#!/usr/bin/env python
#

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import argparse
import logging


### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description='process bed file')
parser.add_argument('-b', '--bed', required=True, 
                   help='bed file to be processed')
parser.add_argument('-o', '--out_bed', required=True, 
                   help='name of output')
##### PARAMs #####
bedFields = ['chrm', 'start', 'stop', 'name', 'score', 'strand']

##### SCRIPT #####
if __name__ == '__main__':
    # load files
    logging.basicConfig(level=logging.INFO)

    args = parser.parse_args()
    logging.info('reading in unprocess bed file: %s'%args.bed)
    with open(args.bed) as f:
        lines = f.readlines()
    
    # go through lines and save fields
    logging.info('processing annotation lines')
    line_dict = {}
    for i, line in enumerate(lines):
        line_series = pd.Series(line.split('\t')[:-1], index=bedFields)
        line_annots = pd.Series({val.split()[0]:val.split()[1].replace('"', '') for val in line.split('\t')[-1].strip('\n;').split('; ')})
        line_dict[i] = pd.concat([line_series, line_annots])
    
    # concat fields
    logging.info('concatenating annotations')
    line_table = pd.concat(line_dict).unstack()
    
    # reorder columns
    cols = line_table.columns.tolist()
    cols_ordered = bedFields + [col for col in cols if col not in bedFields]
    line_table = line_table.loc[:, cols_ordered]

    # save output with fields prepended by comment for use in bed
    logging.info('saving processed output to: %s'%args.out_bed)
    line_table.columns = ['#'+col if i==0 else col for i, col in enumerate(cols_ordered)]
    line_table.to_csv(args.out_bed, sep='\t', index=False)