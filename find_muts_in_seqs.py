#!/usr/bin/env python
#
# i.e.  %run /home/sarah/puflibs/find_muts_in_seqs.py --bed analysis/hPUM2_1muts/clip_signal/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.clip.bed --dat analysis/hPUM2_1muts/sequence/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.fasta.dat

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from fittinglibs import seqfun

# import args
parser = argparse.ArgumentParser()
parser.add_argument('--bed', help='file giving bed file')
parser.add_argument('--dat', help='file giving seq info')
args = parser.parse_args()

# find mutations
bed_fields = ['chrm', 'start', 'stop', 'name',  'score', 'strand',];
bed_data = pd.read_table(args.bed, index_col='name')
bed_data.loc[:, 'seq'] =  pd.read_table(args.dat, header=None,  names=['name', 'seq'], index_col='name', squeeze=True)

# reverse complement if on opposite strand
consensus = 'TGTATATA'
rc_seqs = []
for name, group in bed_data.groupby('strand'):
    if name == '+':
        rc_seqs.append(group.seq.str.upper())
    else:
        rc_seqs.append(pd.Series({key: seqfun.rc(s.upper()) for key, s in group.seq.iteritems()}))
rc_seqs = pd.concat(rc_seqs)

# now for every seq, find muts

