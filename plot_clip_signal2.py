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
dirname = 'analysis/hPUM2_1muts_100kb/clip_signal/plus_strand/'
filename_post = 'hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.100kb.tracks.txt.gz'
filename_dict = {key:os.path.join(dirname, '%s.%s'%(key, filename_post)) for key in ['input', 'rep1', 'rep2']}
data = pd.concat({name:pd.read_csv(filename, compression='gzip', header=None) for name, filename in filename_dict.items()})
dirname = 'analysis/hPUM2_1muts_100kb/clip_signal/minus_strand/'
filename_post = 'hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.100kb.tracks.txt.gz'
filename_dict = {key:os.path.join(dirname, '%s.%s'%(key, filename_post)) for key in ['input', 'rep1', 'rep2']}
data_minus = pd.concat({name:pd.read_csv(filename, compression='gzip', header=None) for name, filename in filename_dict.items()})

# load bed regions
bed_file = 'analysis/hPUM2_1muts_100kb/beds/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.100kb.bed'
bed_fields = ['chrm', 'start', 'stop', 'name',  'score', 'strand',];
bed_data = pd.read_table(bed_file, header=None, names=bed_fields + ['annotation', 'gene', 'exon'], index_col='name')


# find clip enrichment
data_agg2 = (data.sum(axis=1).unstack(level=0).fillna(0) - data_minus.sum(axis=1).unstack(level=0).fillna(0)).astype(float)
clip_enrichment = (data_agg2).apply(lambda x: np.log2((x.rep1 + x.rep2)/(2.*x.input)), axis=1).replace([-np.inf, np.inf], np.nan)
clip_enrichment.index = bed_data.index

# save




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
