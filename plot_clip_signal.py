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

### LOAD files ####
filenames = {'input':'CLIP/hPUM2/consensus_signal//input_plus.hPUM2.0_mut.exons.st.merge.above_0.01_both.st.tracks.txt.gz',
             'rep1':'CLIP/hPUM2/consensus_signal//rep1_plus.hPUM2.0_mut.exons.st.merge.above_0.01_both.st.tracks.txt.gz',
             'rep2':'CLIP/hPUM2/consensus_signal//rep2_plus.hPUM2.0_mut.exons.st.merge.above_0.01_both.st.tracks.txt.gz'}
data = pd.concat({name:pd.read_csv(filename, compression='gzip', header=None) for name, filename in filenames.items()})

# plot integrated signal
data_agg = data.groupby(level=0).mean().transpose()
data_agg.loc[:, 'x'] = np.arange(-100, 101)
data_agg.plot(x='x', figsize=(3.5,3)); plt.ylabel('aggregate signal'); plt.tight_layout()

# order the motif sites
data_agg2 = data.sum(axis=1).unstack(level=0)
data_agg2_log = np.log10(data_agg2)
g = sns.PairGrid(data_agg2_log.loc[:, ['rep1', 'rep2', 'input']]); g.map_upper(plt.scatter, marker='.', edgecolor='none', color='k', s=5);
for ax in g.axes.flat:
    ax.set_xlim([-2, 5]); ax.set_ylim([-2, 5])
order = data_agg2.loc[:, ['rep1', 'rep2']].mean(axis=1).sort_values(ascending=False).dropna().index.tolist()


data_log_zscore = pd.concat([np.log10(mat).apply(lambda x: (x - x.mean())/x.std()) for key, mat in data.groupby(level=0)])
# plot heatmaps
plt.figure(); sns.heatmap(data_log_zscore.loc['rep1'].loc[order])

# compare to RNA seq
peak_id = pd.read_table('RNAseq/motifs/hPUM2/hPUM2.0_mut.exons.st.merge.above_0.01_both.st.peak_id', index_col=0, names=['motif', 'exon', 'gene'], header=None)
peak_id_expanded = {}
for idx, exon, gene in peak_id.itertuples():
    for key in gene.split(','):
        peak_id_expanded[key] = pd.Series([idx, exon], index=['motif', 'exon'])
peak_id_expanded = pd.concat(peak_id_expanded).unstack()

data_agg2_log.index = peak_id.index
data_agg2.index = peak_id.index
data_agg2_log_reindex = pd.DataFrame(index=peak_id_expanded.index, data=data_agg2_log.loc[peak_id_expanded.motif].values, columns=data_agg2_log.columns)
data_agg2_reindex = pd.DataFrame(index=peak_id_expanded.index, data=data_agg2.loc[peak_id_expanded.motif].values, columns=data_agg2.columns)

# load Tpm
tpm = pd.read_table('RNAseq/gene_quant/rna_seq_combined.tpm.above_0.01_both.dat', index_col=0)
log_tpm = np.log10(tpm+1E-3)
clip_enrichment = data_agg2_reindex.apply(lambda x: (x.rep1 + x.rep2)/(2*x.input), axis=1)
