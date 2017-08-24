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
import itertools
import operator
import functools
from hjh import mutations

# import args
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--consensus', help='consensus sequence. Default=TGTATATA', default='TGTATATA')
parser.add_argument('--seq', help='file giving CPannot-type file')
parser.add_argument('--out', help='filename to save to')

args = parser.parse_args()


def make_seq_tag(name, group, attributes):
    """From grouped true/false matrices, make seq tag."""
    cols = [seq for seq, mybool in zip(group.columns.tolist(), name) if mybool]
    if not cols:
        return pd.Series('', index=group.index)
    tags = []
    for seq in cols:
        taglist = [np.power((i+1)*2, df.loc[group.index, seq].astype(int)) for i, df in enumerate(attributes)]
        flags = functools.reduce(operator.mul, taglist, 1)
        tags.append(pd.Series(seq, index=group.index) + '_' + flags.astype(str))
    alltags = pd.Series([';'.join(vec) for vec in zip(*tags)], index=group.index)
    return alltags

def make_seq_tag2(data, attribute_cols=['more_than_one', 'close_to_edge']):
    """From matrix, groupby level 0 """
    alltags = {}
    for idx, group in data.groupby(level=0):
        idx_tags = []
        for (idx2, seq), row in group.iterrows():
            taglist = [np.power((i+1)*2, val) for i, val in enumerate(row.loc[attribute_cols])]
            flag = functools.reduce(operator.mul, taglist, 1)
            idx_tags.append('%s_%d'%(seq, flag))
        alltags[idx] = ';'.join(idx_tags)
    return pd.Series(alltags)

if __name__ == '__main__':
    # load sequence data in chunks. Look for specific sequences
    seqs_to_look = mutations.singles(args.consensus, rna=False) + ['TGTA']
    double_muts = [seq for seq in np.unique(list(itertools.chain(*[mutations.singles(seq, rna=False, includeNoMut=False) for seq in seqs_to_look])))
                   if seq not in seqs_to_look]
    close_threshold = 5               


    chunk = pd.read_table(args.seq, index_col=0)
    # first do for single muts, he  double muts
    find_results = pd.concat({s:chunk.sequence.str.find(s) for s in seqs_to_look}).replace(-1, np.nan).dropna().swaplevel(0, 1).sort_index().astype(int)
    rfind_results = pd.concat({s:chunk.sequence.str.rfind(s) for s in seqs_to_look}).replace(-1, np.nan).dropna().swaplevel(0, 1).sort_index().astype(int)
    
    seq_present = (find_results > -1)
    more_than_one = rfind_results > find_results
    #seq_lengths = chunk.sequence.str.len()
    #right_side_close = pd.Series({(idx, seq):(seq_lengths- close_threshold + len(args.consensus)).loc[idx] for idx, seq  in rfind_results.index.tolist() })
    #close_to_edge = (find_results <= close_threshold)|(rfind_results >= right_side_close)
    
    #data = pd.concat([seq_present.rename('seq_present'), more_than_one.rename('more_than_one'), close_to_edge.rename('close_to_edge')], axis=1)
    data = pd.concat([seq_present.rename('seq_present'), more_than_one.rename('more_than_one')], axis=1)
    tags = data.reset_index(level=1).level_1 + '_' + data.reset_index(level=1).more_than_one.astype(int).astype(str)
    tag_mat = tags.groupby(level=0).apply(lambda x: ';'.join(x))
    #tag_mat = make_seq_tag2(data)
    #tag_mat.sort_index(inplace=True)
    data_all = pd.concat([chunk, tag_mat.rename('seq_tag')], axis=1)
    tags_singles.append(tag_mat.loc[chunk.index])


    
    # save
    if args.out is None:
        args.out = os.path.splitext(args.seq)[0] + '_out' + os.path.splitext(args.seq)[1]
    tags_singles.to_csv(args.out, sep='\t')
     