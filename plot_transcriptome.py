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


# import args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input file')
parser.add_argument('-o', '--output', help='output file')

parser.add_argument('--mode', help='which analysis to run')
args = parser.parse_args()

filenames_all = {}


if args.mode=='plot_single_muts':
    """Find the single muts and plot median dG measured with transcriptome."""
    ddG = pd.read_table('K562_encode/annotations/RNAmap/hPUM2_25.rel_dG.dat', header=None, index_col=0, squeeze=True, names=['name', 'ddG'])
    variant_table = fileio.loadFile('transcriptome/data/transcribable_pum2pos_sub.CPvariant')
    lib_char = fileio.loadFile('transcriptome/final_seqs/transcribable_undet_comb_tag.libChar').loc[variant_table.index]
    ddG.index = ['TGTA_0;' + s + '_0' if s[:4]=='TGTA' else s + '_0' for s in ddG.index]
    
    variant_table.loc[~lib_char.seq_tag.isnull(), 'ddG'] = ddG.loc[lib_char.loc[~lib_char.seq_tag.isnull()].seq_tag].values
    mat = pd.concat([lib_char.annotation, variant_table.dG, variant_table.ddG], axis=1)
    
    y = mat.groupby(['annotation', 'ddG'])['dG'].median()
    yerr = mat.groupby(['annotation', 'ddG'])['dG'].std()/np.sqrt(mat.groupby(['annotation', 'ddG'])['dG'].size())
    data = pd.concat([y.rename('y'), yerr.rename('yerr')], axis=1).reset_index('ddG')
    data.plot(x='ddG', y='y', yerr='yerr', kind='scatter', figsize=(4,4))
    
    plt.xlim([-0.1, 4])
    plt.ylim(-11.5, -7.6)

    # make fasta of consensus seqs
    with open('transcriptome/final_seqs/transcribable_undet_comb_tag_consensus3.fasta', 'w') as f:
        for name, group in lib_char.groupby('seq_tag'):
            if name in [x for x, val in (ddG<0.5).iteritems() if val]:
                for idx, seq in group.sequence.iteritems():
                    f.write('>%s;%d\n'%(name, idx))
                    f.write('%s\n'%seq)
    
    
    
if args.mode == 'remake_fastqs_from_cpseq':
    """Given a cpseq file, remake r1 and r2 fastqs."""
    input_cpseq = args.input
    output_fastq = args.output
    with open(input_cpseq) as f:
        lines = f.readlines()
        
    read1_out = open(output_fastq + '.R1.fastq', 'w')
    read2_out = open(output_fastq + '.R2.fastq', 'w')
    
    for line in lines:
        key, _, read1, read1_qual, read2, read2_qual, _, _, _, _ = line.strip('\n').split('\t')
        read1_out.write('@{key}\n{seq}\n+\n{qual}\n'.format(key=key, seq=read1, qual=read1_qual))
        read2_out.write('@{key}\n{seq}\n+\n{qual}\n'.format(key=key, seq=read2, qual=read2_qual))
        
if args.mode == 'make_cpseq_from_single_fastq':
    """Given a cpseq file, remake r1 and r2 fastqs."""
    input_fastq = args.input
    output_cpseq = args.output
    with open(input_fastq) as f:
        lines = f.readlines()

    with open(output_cpseq, 'w') as f:        
        for key, seq, _, qual in zip(lines[::4], lines[1::4], lines[2::4], lines[3::4]):
            f.write('{key}\t{seq}\t{qual}\t\t\t\t\t\t\t\n'.format(key=key.strip('\n')[1:], seq=seq.strip('\n'), qual=qual.strip('\n')))

    

if args.mode=='find_subset_cpannots':
    """Load CPfitted and variant table and find distributions of dGs"""
    max_num = int(1E4)
    sub_data = []
    for name, group in pd.DataFrame(tags_singles.rename('tags')).groupby('tags'):
        if len(group) > max_num:
            index = np.random.choice(group.index.tolist(), size=max_num, replace=False)
            sub_data.append(group.loc[index])
        else:
            sub_data.append(group)
    sub_data = pd.concat(sub_data)
    sub_data.to_pickle('trimmed_fastqs/Transcribable_S2_L001_bottom.st.sub.CPannot.pkl')
    subset_index = sub_data.index.tolist()
    # find sequences of those in this subset
    sequences = []
    for i, chunk in enumerate(pd.read_table(args.seq, header=None, names=['clusterID', 'sequence'], chunksize=1E6)):
        print i
        sequences.append(chunk.set_index('clusterID').sequence.str.upper())
    sequences = pd.concat(sequences)
    
    sequences_sub = sequences.loc[subset_index]
    
    # for those in this subset, find variant_number/sequence keys
    lib_char = {}
    annotated_clusters = []
    for i, (name, group) in enumerate(pd.DataFrame(sequences_sub).groupby('sequence')):
        annotated_clusters.append(pd.Series(i, index=group.index))
        lib_char[i] = name
    lib_char = pd.DataFrame(pd.Series(lib_char).rename('variant_number'))
    annotated_clusters = pd.DataFrame(pd.concat(annotated_clusters).rename('variant_number'))
    lib_char = pd.concat([lib_char, pd.concat([annotated_clusters, sub_data], axis=1).groupby('variant_number').first()], axis=1)
    
    # for each of the single mutants, find the tag when it is far from the center and there is no more than one
    standard_tags = {}
    for seq in seqs_to_look + ['TGTA']:
        
        if seq == 'TGTA':
            for tgta_flag in [1,2,4,8]:
                tag = 'TGTA_%d'%(tgta_flag)
                standard_tags[tag] = pd.Series([seq, tgta_flag], index=['seq', 'flag'])
        else:
            for seq_append, flag in itertools.product(['TGTA_%d'%i for i in [1,2,4,8]] + [''], [1,2,4,8]):
                if seq_append:
                    if 'TGTA' < seq:
                        tag = '%s;%s_%d'%(seq_append, seq, flag)
                    else:
                        tag = '%s_%d;%s'%(seq, flag, seq_append)
                else:
                    tag = '%s_%d'%(seq, flag)
                standard_tags[tag] = pd.Series([seq, flag], index=['seq', 'flag'])
    standard_tags = pd.concat(standard_tags).unstack()
    
    df = (standard_tags.loc[lib_char.tags])
    df.index = lib_char.index
    lib_char.loc[:, '']
    
    