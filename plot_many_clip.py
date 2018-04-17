#!/usr/bin/env python
#

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import argparse
import itertools
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import functools
import scipy.stats as st
from scikits.bootstrap import bootstrap
from fittinglibs import plotting, seqfun
from tectolibs import tectplots
import scipy.cluster.hierarchy as sch
from puflibs import processing, predictions, seqmodel
# function definitions
def get_log_enrichment_counts(counts1, counts2):
    return np.log2(counts1/counts1.sum() / (counts2/counts2.sum()))


# import args
parser = argparse.ArgumentParser()
parser.add_argument('--mode', help='which analysis to run')
parser.add_argument('--filename_table', help='table of filenames')
args = parser.parse_args()

def find_enrichment(data, num_reads_total=None, method='input'):
    """Determined the fold enrichment by different methods."""
    
    if method=='input':
        weights = 1./num_reads_total*num_reads_total.mean()
        fe = ((data.rep1*weights.rep1 + data.rep2*weights.rep2)/(2*data.input*weights.input)).replace(np.inf, np.nan)
    elif method=='tpm':
        pass
    return re

def process_combined_data(data):
    """Process the data by different methods"""
    data = (data.reset_index().rename(columns={'index':'name'}).groupby(['chrm', 'start', 'stop']).first().reset_index()).copy()
    data = data.groupby(['chrm', 'start', 'stop']).first().reset_index()
    data_expressed = data.loc[data.tpm > 0].copy()
    data_no_neighbors = data.groupby('cluster').applyloc[(data.downstream_bases_to_TGTA > 60)&(data.upstream_bases_toTGTA > 60)].copy()
    data
    
    
if __name__ == '__main__':

    # filenames
    #data = pd.read_table('analysis/output/combined_data.gz', index_col=0, compression='gzip')
    data_0 = pd.read_table('analysis/output/all_unprocessed_st_merged.00.hPUM2_all.random_5e+06.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.temp0.combined_data.01.02.03.04.05.06.07.08.09.ENCFF372VPV.ENCFF141SVY.combined_data.gz', compression='gzip')
    data_37 = pd.read_table('analysis/output/all_unprocessed_st_merged.00.hPUM2_all.random_5e+06.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.temp37.combined_data.01.02.03.04.05.06.07.08.09.ENCFF372VPV.ENCFF141SVY.combined_data.gz',  compression='gzip')
    temperature = 37    
    data = data_37
        
    if args.mode == 'combine_data_all':
        """output of pipeline is split into 10 files. load and combine."""
        filenames = subprocess.check_output('find analysis/output/ -mindepth 2 -name "*combined_data.gz" | sort', shell=True).strip().split()
        data = pd.concat({i:pd.read_table(filename, compression='gzip', index_col=0)
                          for i, filename in enumerate(filenames)})
    
    elif args.mode == 'compare_rna_seq':
        """Make sure the spike in controls of the RNA seq values are quantitative measurements of expression."""
        # note col 'concentration_15' is the concentration of molecules in the version spiked into ENCODE dataset
        nist_ercc_pool15 = pd.read_csv('annotations/nist_ercc_pool15.csv')
        nist_ercc_pool15.index = ['tSpikein_%s'%control for control in nist_ercc_pool15.control]
        rep1 = pd.read_table('analysis/expression/rna_seq_rep1.dat')
        rep2 = pd.read_table('analysis/expression/rna_seq_rep2.dat')
        
        sub_rep1 = rep1.loc[(rep1.transcript_id.str.find('tSpikein_ERCC')==0)].set_index('transcript_id').copy()
        sub_rep1.loc[:, 'actual_concentration'] = nist_ercc_pool15.concentration_15
        
        xlim = np.array([1E-11, 1E-4] )# nmol/ul
        for col in ['expected_count', 'posterior_mean_count', 'TPM', 'FPKM', ]:
            
            min_val = sub_rep1[col].replace(0, np.nan).min()/10.
            plt.figure(figsize=(3,3))
            plt.scatter(sub_rep1.actual_concentration, sub_rep1[col] + min_val, marker='.')
            plt.xscale('log'); plt.yscale('log')
            plt.xlim(xlim)
            plt.ylim(min_val/2, sub_rep1[col].max()*2)
            plt.xlabel('actual concentration (nmol/ul)', fontsize=10)
            plt.ylabel(col, fontsize=10)
            factor = sub_rep1[col].median()/sub_rep1.actual_concentration.median()
            plt.plot(xlim, xlim*factor)
            plt.tight_layout()
        
        # make sure mapping is correct
        tpm_data = pd.concat([rep1.set_index('transcript_id').TPM.rename('rep1'), rep2.set_index('transcript_id').TPM.rename('rep2')], axis=1).reset_index().rename(columns={'transcript_id':'transcript_idx'})
        tpm_data.index = [s.split('.')[0] for s in tpm_data.transcript_idx]
        tpm_combined = np.exp(np.log(tpm_data.loc[:, ['rep1', 'rep2']]).mean(axis=1))
        
        biomart_file = 'annotations/ensemble_gene_converter_biomart.txt'
        biomart_data = pd.read_table(biomart_file, names=['gene_id', 'transcript_id', 'gene_name', 'refseq_id', 'refseq_nc'], header=0)
        biomart_data.loc[:, 'refseq_comb'] = [refseq_id if not str(refseq_id)=='nan' else refseq_nc for idx, refseq_id, refseq_nc in biomart_data.loc[:, ['refseq_id', 'refseq_nc']].itertuples()]

        
        # annotate tpm data with refseq id
        biomart_data.loc[:, 'tpm'] = tpm_combined.loc[biomart_data.transcript_id].values        
        
        
    
    elif args.mode == 'add_col_to_affinity':
        filenames = subprocess.check_output('find analysis/effects/temp_0 -mindepth 2 -name "*affinity.gz" | sort', shell=True).strip().split()
        ddG_no_flip = {}
        ddG_flip = {}
        temperature = 0
        kT = seqmodel.get_ddG_conversion(temperature)
        
        for i, filename in enumerate(filenames):
            seqeffect = pd.read_table(filename, compression='gzip', index_col=0)
            noflip_cols = [idx for idx in seqeffect if idx.find('noflip')==0]
            ddG_no_flip[i] = (kT*np.log(np.exp(seqeffect.loc[:, noflip_cols]/kT).sum(axis=1)))
            flip_cols = [idx for idx in seqeffect if idx not in noflip_cols]
            ddG_flip[i] = (kT*np.log(np.exp(seqeffect.loc[:, flip_cols]/kT).sum(axis=1)))
        ddG_flip = pd.concat(ddG_flip)
        ddG_no_flip = pd.concat(ddG_no_flip)
                          
    elif args.mode == 'compare_window_size':
        """Load the 500 bp window, and calculate enrichment above background for different window sizes"""
        data = pd.read_pickle('analysis/output/all_unprocessed_st_merged.hPUM2.1.random_1e+06.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.combined_data.pkl')

        index_bg = (data.score==0)&(data.tpm > 0.01)
        index_consensus = (data.seq==consensus_seq)&(data.tpm > 0.01)
        index_consensusC = (data.seq==consensus_seqC)&(data.tpm > 0.01)
        index_consensusA = (data.seq==consensus_seqA)&(data.tpm > 0.01)
        consensus_fes = {}
        background_fes = {}
        consensus_counts = {}
        background_counts = {}        
        for interval_size in [10, 20, 40, 80, 160, 320]:

            background_fes[interval_size] = data.loc[index_bg, ['%s_%d'%(key, interval_size) for key in ['rep1', 'rep2']]].sum(axis=1)/data.loc[index_bg].tpm/interval_size
            consensus_fes[interval_size] = data.loc[index_consensus, ['%s_%d'%(key, interval_size) for key in ['rep1', 'rep2']]].sum(axis=1)/data.loc[index_consensus].tpm/interval_size
            background_counts[interval_size] = (data.loc[index_bg, ['%s_%d'%(key, interval_size) for key in ['rep1', 'rep2']]].sum(axis=1)/interval_size).mean()
            consensus_counts[interval_size] = (data.loc[index_consensus, ['%s_%d'%(key, interval_size) for key in ['rep1', 'rep2']]].sum(axis=1)/interval_size).mean()

        consensus_fes = pd.concat(consensus_fes, axis=1)
        consensus_counts = pd.Series(consensus_counts)
        background_counts = pd.Series(background_counts)
        reads_per_bp_80 = consensus_fes.quantile(0.8)
        mean_w_reads = (consensus_fes>0).mean()
        
        g = sns.FacetGrid(pd.concat([reads_per_bp_80.rename('reads_per_bp_80th'), mean_w_reads.rename('fraction_usable')], axis=1).reset_index(), hue='index', palette='viridis');
        g.map(plt.scatter, 'reads_per_bp_80th', 'fraction_usable')

        g = sns.FacetGrid(pd.concat([mean_w_reads.rename('fraction_usable'), (consensus_counts/background_counts).rename('fold_enrichment')], axis=1).reset_index(), hue='index', palette='viridis');
        g.map(plt.scatter, 'fraction_usable', 'fold_enrichment')
        
        clip_signal = (pd.read_csv('analysis/clip/all_unprocessed_st_merged.hPUM2.1.random_1e+06.ann.filt.rep1.ENCFF231WHF.R2.bedGraph.500.tracks.txt.gz', compression='gzip', index_col=0) +
                       pd.read_csv('analysis/clip/all_unprocessed_st_merged.hPUM2.1.random_1e+06.ann.filt.rep2.ENCFF732EQX.R2.bedGraph.500.tracks.txt.gz', compression='gzip', index_col=0))
        xvalues = np.arange(-250, 251)
        num_up_bases = 40
        offset = 15
        plt.figure(); plt.plot(xvalues, clip_signal.loc[index_bg].mean());
        plt.plot(xvalues, clip_signal.loc[index_consensus].mean())
        plt.axvline(-num_up_bases-offset)
        plt.axvline(-num_up_bases-offset+offset+2*num_up_bases)
        plt.errorbar(xvalues, clip_signal.loc[index_consensus].mean(), yerr=clip_signal.loc[index_consensus].std()/np.sqrt(index_consensus.sum())*1.96)
        plt.errorbar(xvalues, clip_signal.loc[index_consensusA].mean(),  yerr=clip_signal.loc[index_consensusA].std()/np.sqrt(index_consensusA.sum())*1.96)
        plt.errorbar(xvalues, clip_signal.loc[index_consensusC].mean(), yerr=clip_signal.loc[index_consensusC].std()/np.sqrt(index_consensusC.sum())*1.96)
    
    elif args.mode == 'apply_model_to_seqs':
        """ """
        basename = 'annotations/RNAmap/qMotif_table_05_012318_1_'
        
        xlim = np.array([0, 5])
        base_params = pd.read_csv(basename + 'term1.csv', index_col=0) #.stack().values
        flip_params = pd.read_csv(basename + 'term2_single.csv', index_col=0) #.stack().values
        dflip_params = pd.read_csv(basename + 'term2_double.csv', index_col=0, squeeze=True) #.stack().values, 10])
        coupling_params= pd.read_csv(basename + 'term3.csv', index_col=0, squeeze=True) #.stack().values
        
        passed_sequence = data.iloc[0].seq_rna
        seqmodel2.additive_PUF_flip_model(passed_sequence, flip_params, base_penalties, coupling_params, double_flip_params, temperature)
    
    
    elif (args.mode == 'plot_flip_annot' or args.mode == 'plot_clip_vs_ddG' or args.mode == 'plot_annotation_vs_ddG' or
          args.mode=='plot_ss_ddG_versus_clip' or args.mode =='plot_versus_proximity'):
        """Using the New model to find effects, find ddG"""
        pass
        data.loc[:, 'flip_annot'] = 'noflip'
        data.loc[data.ddG_noflip - data.ddG > 0.5, 'flip_annot'] = 'flip'
        #data.loc[data.ddG_oneflip - data.ddG > 0.75, 'flip_annot'] = 'dflip'
        
        if args.mode == 'plot_flip_annot':
            """plot the sites that have flip annotation or not."""
            data.loc[:, 'is_random'] = [s.find('hPUM2')!=0 for s in data.name]
            index_subset = np.random.choice((data.loc[(data.tpm>0)&(~data.is_random)].index.tolist()), size=5000, replace=False)
            g = sns.FacetGrid(data=data.loc[index_subset],  hue='flip_annot'); g.map(tectplots.scatter, 'ddG_noflip', 'ddG', marker='.', s=10)
            xlim = np.array([-0.5, 8])
            plt.plot(xlim, xlim, 'k--')
            plt.plot(xlim, xlim-0.5, ':', color='0.5')
            plt.xlim(xlim)
            plt.ylim(xlim)
            data.loc[:, 'dddG_noflip'] = data.ddG_noflip - data.ddG
            plt.figure(figsize=(3,3)); sns.distplot(data.loc[(data.tpm>0)&(~data.is_random)].dddG_noflip, color=sns.color_palette()[0], bins=bins, kde=False)
        
        ddG_binedges = np.hstack([np.linspace(data.ddG.min(), 4.5, 25), data.ddG.max()])
        data.loc[:, 'binned_ddG'] = pd.cut(data.ddG, ddG_binedges,
                                           include_lowest=True, precision=2)
        data.loc[:, 'binned_logtpm'] = pd.cut(np.log10(data.tpm), (np.log10(data.tpm).replace(-np.inf, np.nan).dropna()).quantile(np.linspace(0, 1, 6)), include_lowest=True)
        order = data.groupby('binned_ddG').first().index.tolist()
        
        if args.mode == 'plot_clip_vs_ddG':
            data.loc[:, 'binned_ddG_noflip'] = pd.cut(data.ddG_noflip, ddG_binedges,
                                               include_lowest=True, precision=2)
            med_background_val = data.loc[(data.tpm>0)&(data.ddG > 4.5)].clip_signal_per_tpm.median()
            med_background_val_input = data.loc[(data.tpm>0)&(data.ddG > 4.5)].clip_input_per_tpm.median()
           
            subdata = data.loc[(data.tpm > 0)]
            subdata.loc[:, 'clip_signal_per_tpm_fold'] = subdata.clip_signal_per_tpm/med_background_val
            subdata.loc[:, 'clip_input_per_tpm_fold'] = subdata.clip_input_per_tpm/med_background_val_input

            ylim = [0.5, 100]
            ymax = subdata.groupby('binned_ddG').get_group(order[0])['clip_signal_per_tpm_fold'].median() 
            x_order = pd.Series({name:0.5*(group.ddG.max() + group.ddG.min()) for name, group in subdata.groupby('binned_ddG')})
            for xval, yval in itertools.product(['binned_ddG', 'binned_ddG_noflip'], ['clip_input_per_tpm_fold', 'clip_signal_per_tpm_fold']):

                func = functools.partial(sns.factorplot, data=subdata, hue='flip_annot', y=yval, estimator=np.median,
                                         errwidth=0.5,capsize=1, linestyles='', marker='.', hue_order=['noflip', 'flip'])
                g = func(x=xval); plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
                plt.axhline(1, color='0.5', linestyle='--');
                plt.yscale('log')
                plt.ylim(ylim)
 
                
                y = ymax*np.exp(x_order/seqmodel.get_ddG_conversion(temperature))
                plt.plot(np.arange(len(order)), y.loc[order], 'k--')

                if xval == 'binned_ddG' and yval == 'clip_signal_per_tpm_fold':
                    # add in the ddG_noflip estimations.
                    #noflip_mat = pd.concat({name:pd.cut(group.ddG_noflip, ddG_binedges, precision=2, include_lowest=True).value_counts() for name, group in subdata.loc[subdata.flip_annot=='flip'].groupby('binned_ddG')}).unstack()
                    for name, group in subdata.loc[subdata.flip_annot=='flip'].groupby('binned_ddG'):
                        row = pd.cut(group.ddG_noflip, ddG_binedges, precision=2, include_lowest=True).value_counts()
                        if row.max() > 0:
                            y = group[yval].median()
                            
                            #x = order.index(row.idxmax())
                            #s = float(row.max())/row.sum()*50
                            #plt.scatter(x, y, edgecolor=sns.color_palette()[1], s=s, color='none')
                            for idx_new, val in  (row.astype(float)/row.sum()).iteritems():
                                if val > 0:
                                    x = order.index(idx_new)
                                    s = val*50
                                    plt.scatter(x, y, edgecolor=sns.color_palette()[1], s=s, color='none')
                            
                

            # plot by un flipped dG values           
            x = pd.Series({name:0.5*(group.ddG.max() + group.ddG.min()) for name, group in subdata.groupby('binned_ddG')})
            y = ymax*np.exp(x/seqmodel.get_ddG_conversion(temperature))
            g = func(x='binned_ddG_noflip'); plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
            plt.axhline(med_background_val, color='0.5', linestyle='--')
            plt.plot(np.arange(len(order)), y.loc[order], 'k--')
            ylim = [0.05, 20]
            plt.ylim(ylim)
                
            plt.figure();
            x_noflip = subdata.groupby('flip_annot').get_group('noflip').groupby('binned_ddG')['ddG'].median()
            y_noflip = subdata.groupby('flip_annot').get_group('noflip').groupby('binned_ddG')[yval].median()
            plt.plot(x_noflip, y_noflip, 'o')
            for i, (name, group) in enumerate(subdata.groupby('flip_annot').get_group('flip').groupby('binned_ddG')):
                y = group[yval].median()
                x = group.ddG.median()
                x2 = group.ddG_noflip.median()
                
                color = sns.color_palette()[1]
                plt.plot(x, y, 'o', color=color)
                plt.plot([x, x2], [y, y], '.-', color=color)
            plt.plot(x_noflip, y_noflip[0]*np.exp(x_noflip/seqmodel.get_ddG_conversion(temperature)), 'k--')
            ylim = [0.5, 100]
            plt.ylim(ylim) 
            plt.yscale('log')
            plt.axhline(1, color='0.5', linestyle='--');

                
            

            # plot
            subdata = data.loc[(data.tpm>0)]
            func = functools.partial(sns.factorplot, data=subdata, hue='annotation', hue_order=["5' UTR", "exon", "3' UTR"], y='clip_signal_per_tpm', estimator=np.median,
                                     errwidth=0.5,capsize=1, linestyles='', marker='.')
            g = func(x='binned_ddG'); plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
            plt.axhline(med_background_val, color='0.5', linestyle='--'); plt.yscale('log')
            
            ymax = subdata.groupby('binned_ddG').get_group(order[0]).clip_signal_per_tpm.median()            
            x = pd.Series({name:0.5*(group.ddG.max() + group.ddG.min()) for name, group in subdata.groupby('binned_ddG')})
            y = ymax*np.exp(x/seqmodel.get_ddG_conversion(temperature))
            
            plt.plot(np.arange(len(order)), y.loc[order], 'k--')
            ylim = [0.01, 20]
            plt.ylim(ylim)            
                        
            # plot
            subdata = data.loc[(data.tpm>0)]
            func = functools.partial(sns.factorplot, data=subdata, hue='binned_logtpm',  y='clip_signal_per_tpm', estimator=np.median,
                                     errwidth=0.5,capsize=1, linestyles='', marker='.', palette='viridis')
            g = func(x='binned_ddG'); plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
            plt.axhline(med_background_val, color='0.5', linestyle='--'); plt.yscale('log')

            plt.plot(np.arange(len(order)), y.loc[order], 'k--')
            ylim = [0.01, 20]
            plt.ylim(ylim)
            
            
        # plot by type
        elif args.mode == 'plot_rna_expression_vs_ddG':
            """Compare to RNA expression."""
            subdata = data.loc[(data.tpm>0)&(data.binned_ddG==order[0])].copy()
            subdata.loc[:, 'binned_tpm'] = np.digitize(subdata.tpm, np.logspace(-2, 3, 10))
            pass
        

        elif args.mode == 'plot_annotation_vs_ddG':
            """Using New model to find effects, se how enrichment for UTR changes"""
            expected_fractions = data.loc[data.site.str.find('hPUM') == -1].annotation.value_counts()/(data.site.str.find('hPUM') == -1).sum()
            enrichment = {}
            fraction = {}
            for name, group in data.loc[data.tpm>0].groupby('binned_ddG'):
                
                group1 = group.loc[group.flip_annot == 'noflip']
                group2 = group.loc[group.flip_annot != 'noflip']
                #fracmat[name] = group.annotation.value_counts()/float(len(group))
                fraction[(name, 'noflip')] = (group1.annotation=="3' UTR").mean()
                fraction[(name, 'flip')]   = (group2.annotation=="3' UTR").mean()
                enrichment[(name, 'noflip')] = np.log2((group1.annotation=="3' UTR").mean()/expected_fractions.loc["3' UTR"])
                enrichment[(name, 'flip')]   = np.log2((group2.annotation=="3' UTR").mean()/expected_fractions.loc["3' UTR"])
                
            fraction = pd.Series(fraction).unstack().loc[order]
            fraction.plot(marker='o', linestyle='none', figsize=(3,3))
            plt.xticks(np.arange(len(order)), order, rotation=90)
            plt.axhline(expected_fractions.loc["3' UTR"], linestyle='--', color='k')
            plt.xlim(-1, len(order))
            plt.ylim(0, 1)
            
        elif args.mode == 'plot_versus_proximity':
            """Plot how the distribution of clip enrichment changes if there are nearby UGUA sites"""
            data.loc[:, 'min_distance_to_TGTA'] = data.loc[:, ['upstream_bases_to_TGTA', 'downstream_bases_to_TGTA']].min(axis=1)
            binedges = [0,  60,  150, 251]
            data.loc[:, 'binned_min_dist'] = pd.cut(data.min_distance_to_TGTA, binedges, include_lowest=True)
            
            subdata = data.loc[(data.tpm>0)&(data.annotation=="3' UTR")]
            g = sns.factorplot(data=subdata, x='binned_ddG', hue='binned_min_dist', y='clip_input_per_tpm', estimator=np.median,
                                     errwidth=1,capsize=1, linestyles='', marker='.', palette='viridis')
            plt.yscale('log'); plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
            ymax = subdata.groupby('binned_ddG').get_group(order[0]).clip_signal_per_tpm.median()            
            x = pd.Series({name:0.5*(group.ddG.max() + group.ddG.min()) for name, group in subdata.groupby('binned_ddG')})
            y = ymax*np.exp(x/seqmodel.get_ddG_conversion(temperature))
            
            plt.plot(np.arange(len(order)), y.loc[order], 'k--')
            ylim = [0.01, 20]
            plt.ylim(ylim)
            
        elif args.mode == 'plot_ss_ddG_versus_clip':
            
            subdata = data.loc[(data.binned_ddG==order[0])&(data.tpm>0)]
            subdata.loc[:, 'binned_ddG_ss'] = pd.cut(data.ss_ddG,
                                               #data.ddG.quantile(np.linspace(0, 1, 100)),
                                               np.hstack([np.linspace(subdata.ss_ddG.min(), 20, 25),
                                                         subdata.ss_ddG.max()]),
                                               include_lowest=True, precision=2)
            order_ss = subdata.groupby('binned_ddG_ss').first().index.tolist()
            g = sns.factorplot(data=subdata, x='binned_ddG_ss', y='clip_signal_per_tpm', estimator=np.median,
                                     errwidth=1,capsize=2, linestyles='', marker='.', palette=['k'])
            plt.yscale('log'); plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)            
            ylim = plt.gca().get_ylim()
            
            # plot line
            ymax = subdata.groupby('binned_ddG_ss').get_group(order_ss[0]).clip_signal_per_tpm.median()            
            x = np.linspace(subdata.ss_ddG.min(), 20)
            y = ymax*np.exp(x/seqmodel.get_ddG_conversion(temperature))
            
            plt.plot(x, y, 'k--')
            plt.ylim(ylim)
    
    elif args.mode == 'plot_indiv_sites':
        """Rather than aggregate behavior, examine individual sites"""
        norad_refseqid = 'NR_027451'
        
        counts = {}
        for i in np.arange(11):
            print i
            directory = 'analysis/clip/split_%d/'%i
            filenames =  [os.path.join(directory, filename) for filename in os.listdir(directory) if filename[-6:]=="txt.gz"]
            keys = [filename.split('.')[-8] for filename in filenames]
            for key, filename in zip(keys, filenames):
                data_table = pd.read_csv(filename, compression='gzip', index_col=0)
                counts[(i, '%s'%(key))] = processing.get_counts_from_counts_table(data_table, )
        counts = pd.concat(counts).unstack(level=1)        
        clip_signal_counts = counts.rep1 + counts.rep2
        
        data = pd.concat([data.set_index(['split_id', 'name']), clip_signal_counts.rename('clip_signal_counts')], axis=1).reset_index()

        norad_subset = data.loc[data.refseq_id==norad_refseqid].sort_values('start').copy()
        dist_to_next_site = pd.Series({idx:next_start - stop  for idx, stop, next_start in zip(norad_subset.index.tolist()[:-1], norad_subset.stop.iloc[:-1], norad_subset.start.iloc[1:])})
        
        # for all sites within 40 bp, take the max
        for idx, dist in dist_to_next_site.iteritems():
            pass
            
            
    elif args.mode == 'compare_to_pum12_kd':
        """Load supp data from NAR paper and compare sites."""
        
        expression_data = pd.read_csv('annotations/nar-01280/supp_table4.csv')
        expression_data.loc[:, 'lfc'] = expression_data.lfc.replace('#NAME?', np.nan).astype(float).replace(np.inf, np.nan)        
        expression_data.loc[:, 'log_fpkm_cntrl'] = np.log10(expression_data.FPKM_NTC)

        # choose randomly from a subset of similarly expressed genes in the control
        row_subsets = {}
        length_factor = 5
        for name, subset in zip(['up_matched', 'down_matched'], [expression_data.sig_up, expression_data.sig_down]):

            kernel = st.gaussian_kde(expression_data.loc[subset].log_fpkm_cntrl)
            prob_density = kernel(expression_data.loc[~subset].log_fpkm_cntrl)

            kernel_start = st.gaussian_kde(np.random.choice(expression_data.loc[~subset].log_fpkm_cntrl.replace([-np.inf], np.nan).dropna(), 1000))
            prob_density_orig = kernel_start(expression_data.loc[~subset].log_fpkm_cntrl)

            weights = (pd.Series(prob_density/prob_density_orig/
                                np.nansum((prob_density/prob_density_orig)), index=expression_data.loc[~subset].index).
                       replace([-np.inf, np.inf], np.nan).fillna(0))

            row_subsets[name] = np.random.choice(expression_data.loc[~subset].index.tolist(),
                                                  p=weights,
                                                  size=int(subset.sum()*length_factor),
                                                  replace=False)
            row_subsets[name.split('_')[0]] = expression_data.loc[subset].index.tolist()
        
        
        # load the biomart ref
        biomart_data = pd.read_table('annotations/ensemble_gene_converter_biomart.txt', names=['gene_id', 'transcript_id', 'gene_name', 'refseq_id', 'refseq_nc'], header=0)
        biomart_data.loc[:, 'refseq_comb'] = [refseq_id if not str(refseq_id)=='nan' else refseq_nc for idx, refseq_id, refseq_nc in biomart_data.loc[:, ['refseq_id', 'refseq_nc']].itertuples()]
        data.loc[:, 'gene_name'] = pd.Series(biomart_data.groupby('refseq_comb').first().loc[data.refseq_id.dropna()].gene_name.values, index=data.refseq_id.dropna().index)
        

        
        # group by the gene and find occupancy and other metrics
        RT = -seqmodel.get_ddG_conversion(temperature=37)
        ss_ddG_threshold = 10 # kcal/mol
        occupancy_data = {}
        for name, group in data.groupby('gene_name'):
            # filter all for ss structure
            group_3UTR = group.loc[(group.annotation=="3' UTR")&(group.ss_ddG < ss_ddG_threshold)]
            occupancy_3UTR = np.exp(-group_3UTR.ddG/RT).sum()
            occupancy_noflip_3UTR = np.exp(-group_3UTR.ddG_noflip/RT).sum()
            num_consensus_3UTR = (group_3UTR.ddG < 0.5).sum()
            num_sites_2kc_3UTR = (group_3UTR.ddG < 2).sum()
            num_sites_between1and4kc_3UTR = ((group_3UTR.ddG < 4)&(group_3UTR.ddG >= 1)).sum()
            min_dG_3UTR = group_3UTR.ddG.min()
            occupancy_not3UTR = np.exp(-group.loc[group.annotation!="3' UTR"].ddG/RT).sum()
            occupancy_data[name] = pd.Series({'occupancy_3UTR':occupancy_3UTR,
                                              'occupancy_noflip_3UTR':occupancy_noflip_3UTR,
                                              'occupancy_not3UTR':occupancy_not3UTR,
                                              'min_ddG':min_dG_3UTR,
                                              'num_consensus_3UTR':num_consensus_3UTR,
                                              'num_sites_2kc_3UTR':num_sites_2kc_3UTR,
                                              'num_sites_between1and4kc_3UTR':num_sites_between1and4kc_3UTR})
        occupancy_data = pd.concat(occupancy_data).unstack()
        occupancy_data.loc[occupancy_data.min_ddG.isnull(), 'min_ddG'] = 10 #kcal.mol
        # find roc
        roc_curves = {}
        for name in occupancy_data:
            
            expression_data.loc[:, 'occupancy'] = occupancy_data.loc[expression_data.gene, name].values
            expression_data_sub = expression_data.dropna(subset=['lfc', 'occupancy']).copy()
            if name=='min_ddG':
                predicted_up=False
            else:
                predicted_up = True
            #roc_curves[(name)] =  processing.find_roc_data(expression_data_sub, 'occupancy', 'sig_up')
            for expression_threshold in [-np.inf, 1, 1.5 ]:
                roc_curves[(expression_threshold, name)] =  processing.find_roc_data(expression_data_sub.loc[expression_data_sub.log_fpkm_cntrl >= expression_threshold], 'occupancy', 'sig_up', predicted_up=predicted_up, )
        roc_curves = pd.concat(roc_curves, names=['cntrl_expression', 'occ_def', 'occ_val'])
        g = sns.FacetGrid(data=roc_curves.reset_index(), hue='occ_def', col='cntrl_expression', hue_order=['occupancy_3UTR', 'occupancy_noflip_3UTR', 'occupancy_not3UTR'], palette=['b', 'r', '0.5'], ); g.map(plt.plot, 'fpr', 'tpr', )
        g = sns.FacetGrid(data=roc_curves.reset_index(), hue='occ_def', col='cntrl_expression', hue_order=['occupancy_3UTR', 'min_ddG', 'num_sites_2kc_3UTR', 'num_sites_between1and4kc_3UTR'], palette=['b', 'g', 'm', 'c'], ); g.map(plt.plot, 'fpr', 'tpr', )

        # find ROC curve for genes with NO consensus sites
        expression_data.loc[:, 'occupancy'] = occupancy_data.loc[expression_data.gene, 'occupancy_3UTR'].values
        expression_data.loc[:, 'has_consensus'] = occupancy_data.loc[expression_data.gene, 'num_consensus_3UTR'].values > 0
        expression_data_sub = expression_data.dropna(subset=['lfc', 'occupancy'])
        roc_noconsesus = processing.find_roc_data(expression_data_sub.groupby('has_consensus').get_group(False), 'occupancy', 'sig_up')
        roc_wconsesus = processing.find_roc_data(expression_data_sub.groupby('has_consensus').get_group(True), 'occupancy', 'sig_up')
        
        g = sns.FacetGrid(roc_curves.loc[-np.inf].reset_index(), hue='occ_def', hue_order=['occupancy_3UTR', 'occupancy_not3UTR'], palette=['b', '0.5']); g.map(plt.plot, 'fpr', 'tpr')
        plt.plot(roc_noconsesus.fpr, roc_noconsesus.tpr, 'b--')
        idx_noconsensus = np.abs(roc_noconsesus.reset_index().loc[:, 'index'] - 1).sort_values().index[0]
        plt.scatter(roc_noconsesus.reset_index().loc[idx_noconsensus].fpr, roc_noconsesus.reset_index().loc[idx_noconsensus].tpr, c='b')
        #plt.plot(roc_wconsesus.fpr, roc_wconsesus.tpr, 'b:')
        idx_gene_body = np.abs(roc_curves.loc[-np.inf].loc['occupancy_not3UTR'].reset_index().occ_val-1).sort_values().index[0]
        plt.scatter(roc_curves.loc[-np.inf].loc['occupancy_not3UTR'].reset_index().loc[idx_gene_body].fpr, roc_curves.loc[-np.inf].loc['occupancy_not3UTR'].reset_index().loc[idx_gene_body].tpr, c='0.5')
        plt.plot([0, 1], [0, 1], 'k--')
        
        # plot the roc plot for the other model
        expression_data.loc[:, 'occupancy'] = occupancy_data.loc[expression_data.gene, 'occupancy_noflip_3UTR'].values
        expression_data_sub = expression_data.dropna(subset=['lfc', 'occupancy'])
        roc_noconsesus_noflip = processing.find_roc_data(expression_data_sub.groupby('has_consensus').get_group(False), 'occupancy', 'sig_up')
        plt.plot(roc_noconsesus_noflip.fpr, roc_noconsesus_noflip.tpr, 'r--')
        

        print processing.get_tpr_fpr(expression_data_sub.has_consensus, expression_data_sub.sig_up)
        print processing.get_tpr_fpr(expression_data_sub.loc[~expression_data_sub.has_consensus].occupancy >= 1, expression_data_sub.loc[~expression_data_sub.has_consensus].sig_up)
        
        roc = processing.find_roc_data(expression_data_sub, 'occupancy', 'sig_up')
        
        # laod 3' UTR length from refseq database
        
        bed_data = pd.read_table('annotations/refseq/hg38_refGene.3UTR.bed', names=variables.bed_fields + ['gene_name'], header=0)
        bed_data.loc[:, 'utr_length'] = bed_data.stop - bed_data.start
        
        expression_data.loc[:, "utr_length_min"] = bed_data.groupby('gene_name')['utr_length'].min().loc[expression_data.gene].values
        expression_data.loc[:, "utr_length_max"] = bed_data.groupby('gene_name')['utr_length'].max().loc[expression_data.gene].values
        expression_data.loc[:, "utr_length_median"] = bed_data.groupby('gene_name')['utr_length'].median().loc[expression_data.gene].values
        roc_curves_length = {}
        for name in ['utr_length_min', 'utr_length_max', 'utr_length_median']:
            roc_curves_length[name] = processing.find_roc_data(expression_data.dropna(subset=['lfc', name]).copy(), name, 'sig_up')
        roc_curves_length = pd.concat(roc_curves_length, names=['length_def', 'length_val'])
        sys.exit()
        
        # for each gene in the expression data, find the sites that are associated with it
        data_subsets = {}
        for key, idxs in row_subsets.items():
            print key
            for i, (idx, gene) in enumerate(expression_data.loc[idxs].gene.iteritems()):
                if i%10 == 0:
                    print i
                possible_refseq_ids = biomart_data.loc[(biomart_data.gene_name==gene)].refseq_comb.dropna().unique().tolist()
                data_subset = data.loc[pd.Series(np.in1d(data.refseq_id, possible_refseq_ids), index=data.index)&
                                       (data.annotation=="3' UTR")].copy()
                
                data_to_save = pd.Series({'min_ddG':data_subset.ddG.min(), 'median_ddG':data_subset.ddG.median(),
                                          'num_sites':len(data_subset), 'num_sites_2kc':(data_subset.ddG < 2).sum(),
                                        'sites_index':data_subset.index.tolist()})
                data_subsets[(key, gene)] = data_to_save
        
    if args.mode == "compare_to_clip_peaks":
        """Compare the motifs found to the clip peaks analysis"""
        clip_peaks = processing.load_bed('analysis/beds/ENCFF372VPV.rep1.ENCFF141SVY.rep2.bed')
        gc_content = pd.read_table('analysis/beds/ENCFF372VPV.rep1.ENCFF141SVY.rep2.bed.nuc', usecols=['4_usercol', '8_pct_gc'], index_col='4_usercol', squeeze=True)
        peak_count = pd.read_table('analysis/beds/ENCFF372VPV.rep1.ENCFF141SVY.rep2.bed.bam_counts', index_col=0)
        total_counts = {'rep2':np.loadtxt('analysis/bams/rep2.ENCFF732EQX.R2.bamcount')}
        which_motifs = pd.read_table('analysis/beds/ENCFF372VPV.rep1.ENCFF141SVY.rep2.bed.peak_id', skiprows=1, header=None, names=['clip_peak', 'motif_site']).replace('.', np.nan)
        peak_score = (peak_count.rep1 + peak_count.rep2)/(peak_count.input + 1)
        num_sites = (which_motifs.dropna().groupby('clip_peak').size()).loc[clip_peaks.name].fillna(0)
        
        site_info = data.loc[data.split_id==10].set_index('name')
        
    ####### OLD 2/2/18 #######
    
    if args.mode == 'compare_consensus_ss_energies':
        """Make plots for the comparing the clip enrichment to the ss energy difference."""
        
        order = ['ss_50', 'ss_75', 'ss_100']
        
        data = pd.concat([motif_data.loc[(motif_data.score>11), ['fold_enrichment']+order]], axis=1).dropna()
        data_melt = pd.melt(data, id_vars=['fold_enrichment'])
        xlim = np.array([0, 20])
        g = sns.FacetGrid(data=data_melt, col='fold_enrichment', xlim=xlim, palette=['0.5']); g.map(plt.scatter, 'value', 'fold_enrichment', marker='.', edgecolor='k')
    
        for i, ax in enumerate(g.axes.flat):
            slope, intercept, rvalue, pvalue, stderr = st.linregress(data.loc[:, order[i]], data.loc[:, 'clip_enrichment'])
            ax.plot(xlim, xlim*slope+intercept, c='k')
            ax.annotate('rsq = %4.3f'%rvalue**2, xy=(0.95, 0.95), xycoords='axes fraction',
                                    horizontalalignment='right',
                verticalalignment='top',)
    
        # plot mean.
        ddG_binned = pd.Series(np.digitize(motif_data.loc[:, order].mean(axis=1), [2, 4, 6, 8, 10]), index=motif_data.index)
        motif_data.loc[:, 'ss_binned'] = ddG_binned
        motif_data.loc[:, 'ss_ddG'] = ddG_binned
    
        # plot
        subset_index = motif_data.ddG < 0.5
        sub_data = motif_data.loc[subset_index].copy()
        ddG_ss = sub_data.loc[subset_index].groupby('ss_binned')[order].median().mean(axis=1)
        
        # take subset with low ddG (affinity)
        cols = ['input', 'rep1', 'rep2']
    
        # do it differently
        num_reads = random_data.loc[:, ['input', 'rep1', 'rep2']].mean()
        #num_reads = pd.Series(1, index=cols)
    
        counts = sub_data.groupby('ss_binned')[cols].mean()/num_reads
        counts_bounds = pd.concat({(name, col):pd.Series(np.abs(bootstrap.ci(group.loc[:, col]/num_reads.loc[col], np.mean) - counts.loc[name, col]),
                                                         index=['eminus', 'eplus'])
                         for col in cols for name, group in sub_data.groupby('ss_binned')})
        counts_err = pd.concat({name:group.unstack().loc[:, name].unstack()
                                for name, group in counts_bounds.groupby(level=2)})
        mat = pd.concat([ddG_ss.rename('ddG'), counts, ], axis=1)
        
        # plot
        plt.figure(figsize=(4,3))
        for col in cols:
            plt.scatter(mat.ddG, mat.loc[:, col], label=col)
            plt.errorbar(mat.ddG, mat.loc[:, col], yerr=[counts_err.loc[key].loc[mat.index.tolist(), col] for key in ['eminus', 'eplus']], fmt=',')
        plt.yscale('log')
        plt.xlabel('$\Delta \Delta G$ (kcal/mol)')
        plt.ylabel('num. reads per site')
        plt.tight_layout()
        
        # plot the 
        g = sns.FacetGrid(sub_data_melt, hue='ss_binned', col='variable'); g.map(tectplots.plot_cdf, 'value'); plt.xscale('log')
    
    
    if args.mode=='find_num_in_annotations':
        """Find the enrichmnet for various annotations"""
        order = ["5' UTR", "exon", "3' UTR"]
        sub_data = motif_data.loc[ motif_data.ddG < 0.5].copy()
        fraction_consensus = sub_data.annotation.value_counts()/len(sub_data)
    
        expected_fraction  = random_data.annotation.value_counts()/len(random_data)
        
        # plot fraction
        tectplots.figure(); sub_data.annotation.value_counts().loc[order].plot(kind='bar'); plt.ylabel('number of sites'); plt.tight_layout()
        
        # plot enrichment
        tectplots.figure(); np.log2(fraction_consensus/expected_fraction).loc[order].plot(kind='bar'); plt.ylabel('log2 obs/expected');plt.tight_layout()
        
        # plot signal in each annotation
        sub_data_melt = pd.melt(sub_data.loc[:, ['annotation'] + cols], id_vars=['annotation'])
        g = sns.FacetGrid(sub_data_melt, hue='annotation', col='variable'); g.map(tectplots.plot_cdf, 'value'); plt.xscale('log')
        g = sns.FacetGrid(sub_data, hue='annotation'); g.map(tectplots.plot_cdf, 'fold_enrichment'); plt.xscale('log')
        
        
        
    if args.mode=='process_meadured_dGs':
        """Load the measured_dG tables and average"""
        usecols = ['dG', 'PUF_KERNEL_SEQ', 'sequence', 'foldFreeEnergy']
        filenames = ['annotations/RNAmap/hPUM2_AM4NT_25.csv', 'annotations/RNAmap/hPUM2_AMEMH_25.csv']
        data = (pd.concat({i:pd.read_csv(filename, usecols=usecols, index_col='sequence')
                          for i, filename in enumerate(filenames)}, names=['rep']).
                rename(columns={'PUF_KERNEL_SEQ':'motif',  'foldFreeEnergy':'dG_fold'}))
        data.loc[:, 'dG_both'] = data.dG + data.dG_fold/3.
        consensus = 'TGTATATA'
    
        ddG = data.groupby('motif')['dG_both'].mean() - data.groupby('motif')['dG_both'].mean().loc[consensus]
        #save: 'annotations/RNAmap/hPUM2_25.rel_dG.dat'
        """
        index_list = []
        for s in ddG.index.tolist():
            vec = [(i, s2) for i, (s1, s2) in enumerate(zip(consensus, s)) if s1!=s2]
            if len(vec)==1:
                index_list.append(vec[0])
            elif len(vec) == 0:
                index_list.append((-1, ''))
            else:
                print 'double mutant found!'
                break
        ddG.index = pd.MultiIndex.from_tuples(index_list)
        """
    
    if args.mode=='compare_to_dGs_measured' or args.mode=='compare_to_ss_dGs' or args.mode=='compare_transcripts_with_binding_sites' or args.mode=='compare_transcripts_more':
        """Compare the measured dGs to the clip enrichmnet."""
        bed_data = motif_data.copy()
    
        
        ddG = pd.read_table('annotations/RNAmap/hPUM2_25.rel_dG.dat', header=None, index_col=0, squeeze=True, names=['name', 'ddG'])
        bed_data.loc[:, 'ddG'] = ddG.loc[bed_data.seq].values
        bed_data.loc[:, 'binned_ddG'] = np.digitize(bed_data.ddG, [0.5, 1.1, 1.7, 2.3])
        bed_data.loc[:, 'log_total_counts'] = np.log10(bed_data.rep1 + bed_data.rep2).replace(-np.inf, np.nan)
        bed_data.loc[:, 'log_tpm'] = np.log10(bed_data.tpm).replace(-np.inf, np.nan)
    
        random_data.loc[:, 'log_total_counts'] = np.log10(random_data.rep1 + random_data.rep2).replace(-np.inf, np.nan)
        random_data.loc[:, 'log_tpm'] = np.log10(random_data.tpm).replace(-np.inf, np.nan)
        
        data_flip = pd.read_pickle('analysis/output/hPUM2_flip_5U.0.input.ENCFF786ZZB.R2.rep2.ENCFF732EQX.rep1.ENCFF231WHF.combined_data.pkl')
        bed_data_flip = data_flip.loc['motif'].copy()
        random_data_flip = data_flip.loc['random'].copy()
        bed_data_flip.loc[:, 'log_total_counts'] = bed_data_flip.apply(lambda x: np.log10(x.rep1 + x.rep2), axis=1).replace(-np.inf, np.nan)
        bed_data_flip.loc[:, 'log_tpm'] = bed_data_flip.apply(lambda x: np.log10(x.tpm), axis=1).replace(-np.inf, np.nan)
        random_data_flip.loc[:, 'log_total_counts'] = random_data_flip.apply(lambda x: np.log10(x.rep1 + x.rep2), axis=1).replace(-np.inf, np.nan)
        random_data_flip.loc[:, 'log_tpm'] =  random_data_flip.apply(lambda x: np.log10(x.tpm), axis=1).replace(-np.inf, np.nan)
            
        if args.mode == 'compare_to_dGs_measured':
            
            # plot log counts
            num_in_bin = bed_data.groupby('ddG').size()
            bed_data_red = pd.concat({'value':bed_data.groupby('ddG')['log_total_counts'].mean(),
                                      'value_err':bed_data.groupby('ddG')['log_total_counts'].std()/np.sqrt(num_in_bin)}, axis=1)
            #bed_data_red.reset_index().plot(x='ddG', y='value', yerr='value_err', marker='o')
            
            # plot tpm avg
            bed_data_red.loc[:, 'tpm_avg'] = bed_data.groupby('ddG')['log_tpm'].mean()
            bed_data_red.loc[:, 'tpm_avg_err'] = bed_data.groupby('ddG')['log_tpm'].std()/np.sqrt(num_in_bin)
            
            # subtract off tpm
            bed_data_red.loc[:, 'value_norm'] =  bed_data_red.value - bed_data_red.tpm_avg + bed_data_red.tpm_avg.mean()
            bed_data_red.loc[:, 'value_norm_err'] = np.sqrt(bed_data_red.value_err**2 + bed_data_red.tpm_avg_err**2)
            
            # find background
            background_value = random_data.log_total_counts.mean()
            background_value_err = random_data.log_total_counts.std()/np.sqrt(len(random_data.log_total_counts.dropna()))
            
            # subtract off background
            bed_data_red.loc[:, 'value_sub'] = bed_data_red.value_norm - background_value
            bed_data_red.reset_index().plot(x='ddG', y='value_sub', yerr='value_norm_err', marker='o', figsize=(3,3), color='r')
                    
            plt.axhline(0, color='k')
            plt.axhline(0-background_value_err, color='k', linestyle=':')
            plt.axhline(0+background_value_err, color='k', linestyle=':')
            xlim = np.array([-0.1, 4])
            plt.xlim(xlim)
            
            ticklocs = np.log10([1, 2, 5, 10])
            plt.yticks(ticklocs, [1, 2, 5, 10])
            
            #plt.plot(xlim, -xlim/np.log(10)/RT + bed_data_red.value_norm.max() , color='k', linestyle='--' )
            plt.ylim(-0.1, 1)
            
            # plot the flipped data
            value = bed_data_flip['log_total_counts'].mean()
            value_err = bed_data_flip['log_total_counts'].std()/np.sqrt(len(bed_data_flip))
    
            background_value = random_data_flip.log_total_counts.mean()
            #background_value_err = random_data_flip.log_total_counts.std()/np.sqrt(len(random_data_flip))
            value_sub = value - background_value
    
            plt.axhline(value_sub, color='r')
            plt.axhline(value_sub+value_err, color='r', linestyle=':')
            plt.axhline(value_sub-value_err, color='r', linestyle=':')
            # find 
        if args.mode == 'compare_to_ss_dGs':
            binedges = bed_data.ss_75.quantile(np.linspace(0, 1, 12))
            binedges.iloc[-1] += 1
            
            binedges = [2, 4, 6, 8, 10]
            binedges = np.arange(1, 12)
            bed_data.loc[:, 'ss_ddG_binned'] = np.digitize(bed_data.ss_75, binedges)
            
            bed_data.loc[:, 'ss_ddG'] = bed_data.groupby('ss_ddG_binned')['ss_75'].mean().loc[bed_data.ss_ddG_binned].values
            bed_data_sub = bed_data.loc[bed_data.ddG < 0.5]
            # plot log counts
            num_in_bin = bed_data_sub.groupby('ss_ddG').size()
            bed_data_red = pd.concat({'value':bed_data_sub.groupby('ss_ddG')['log_total_counts'].mean(),
                                      'value_err':bed_data_sub.groupby('ss_ddG')['log_total_counts'].std()/np.sqrt(num_in_bin)}, axis=1)
            #bed_data_red.reset_index().plot(x='ddG', y='value', yerr='value_err', marker='o')
            
            # plot tpm avg
            bed_data_red.loc[:, 'tpm_avg'] = bed_data_sub.groupby('ss_ddG')['log_tpm'].mean()
            bed_data_red.loc[:, 'tpm_avg_err'] = bed_data_sub.groupby('ss_ddG')['log_tpm'].std()/np.sqrt(num_in_bin)
            
            # subtract off tpm
            bed_data_red.loc[:, 'value_norm'] =  bed_data_red.value - bed_data_red.tpm_avg + bed_data_red.tpm_avg.mean()
            bed_data_red.loc[:, 'value_norm_err'] = np.sqrt(bed_data_red.value_err**2 + bed_data_red.tpm_avg_err**2)
            
            # find background
            background_value = random_data.log_total_counts.mean()
            background_value_err = random_data.log_total_counts.std()/np.sqrt(len(random_data.log_total_counts.dropna()))
    
            # subtract off background
            bed_data_red.loc[:, 'value_sub'] = bed_data_red.value_norm - background_value
            bed_data_red.reset_index().plot(x='ss_ddG', y='value_sub', yerr='value_norm_err', marker='o', figsize=(3,3), color='r')
            
            
            plt.axhline(0, color='k')
            plt.axhline(-background_value_err, color='k', linestyle=':')
            plt.axhline(background_value_err, color='k', linestyle=':')
            xlim = np.array([-0.1, 14])
            plt.xlim(xlim)
            
            ticklocs = np.log10([1, 2, 5, 10])
            plt.yticks(ticklocs, [1, 2, 5, 10])
            plt.ylim(-0.1, 1)
            
            bed_data_red.reset_index().plot(x='ss_ddG', y='value_sub', yerr='value_norm_err', marker='o', figsize=(3,3), color='r')
            ticklocs = np.log10([2.5, 5, 10, 20]) - background_value
            plt.yticks(ticklocs, [2.5, 5, 10, 20])
            plt.ylim(-0.1, 1)
            
            #plt.plot(xlim, -xlim/np.log(10)/RT + bed_data_red.value_norm.max() , color='k', linestyle='--' )
    
            
        if args.mode == 'compare_to_dGs_measured_old':
    
            # bin the ddGs and plot as a violinplot
            sns.factorplot(data=bed_data, x='binned_ddG', y='fold_enrichment', kind='box', notch=True, palette=['0.5'], aspect=0.7); plt.yscale('log'); plt.ylim(0.1, 1000)
            sns.factorplot(data=random_data, x='score', y='fold_enrichment', kind='box', notch=True, palette=['0.5'], aspect=0.7); plt.yscale('log'); plt.ylim(0.1, 1000)
    
            tpm_sum = motif_data.groupby('seq')['tpm'].sum()
            cols = ['input', 'rep1', 'rep2']
        
            # do it differently
            #num_reads = random_data.loc[:, ['input', 'rep1', 'rep2']].mean()
            #num_reads = pd.Series(1, index=cols)
    
            counts = motif_data.groupby('seq')[cols].mean()*weights
            counts_bounds = pd.concat({(name, col):pd.Series(np.abs(bootstrap.ci(group.loc[:, col]*weights.loc[col], np.mean) - counts.loc[name, col]),
                                                             index=['eminus', 'eplus'])
                             for col in cols for name, group in motif_data.groupby('seq')})
            counts_err = pd.concat({name:group.unstack().loc[:, name].unstack()
                                    for name, group in counts_bounds.groupby(level=2)})
            mat = pd.concat([ddG, counts, tpm_sum], axis=1)
            
            # plot
            plt.figure(figsize=(4,3))
            for col in cols:
                plt.scatter(mat.ddG, mat.loc[:, col], label=col)
                plt.errorbar(mat.ddG, mat.loc[:, col], yerr=[counts_err.loc[key].loc[mat.index.tolist(), col] for key in ['eminus', 'eplus']], fmt=',')
            plt.yscale('log')
            plt.xlabel('$\Delta \Delta G$ (kcal/mol)')
            plt.ylabel('num. reads per site')
            plt.tight_layout()
            plt.ylim(1, 500)
    
            plt.figure(figsize=(4,3))
            for col in cols:
                plt.scatter(mat.ddG, mat.loc[:, col], label=col)
                plt.errorbar(mat.ddG, mat.loc[:, col], yerr=[counts_err.loc[key].loc[mat.index.tolist(), col] for key in ['eminus', 'eplus']], fmt=',')
            plt.yscale('log')
            plt.xlabel('$\Delta \Delta G$ (kcal/mol)')
            plt.ylabel('num. reads per site')
            plt.tight_layout()
            plt.ylim(1, 500)
            
            # plot tpm
            mat_melt = pd.melt(mat, id_vars=['ddG', 'tpm'])
            g = sns.FacetGrid(mat_melt, col='variable'); g.map(plt.scatter, 'tpm', 'value')
            
            """
            # plot the enrichment over input
            plt.figure(figsize=(4,3))
            for col in ['rep1', 'rep2']:
                a = mat.loc[:, col]
                b = mat.input
                y = a/b
                yerr = []
                for key in ['eminus', 'eplus']:
                    da = counts_err.loc[key, col]
                    db = counts_err.loc[key].input
                    yerr.append(np.sqrt(((da/a)**2 + (db/b)**2).astype(float))*y)
    
                plt.scatter(mat.ddG, y)
                plt.errorbar(mat.ddG, y, yerr=yerr, fmt=',')
    
            plt.yscale('log')
            plt.ylim(0.5, 40)
            plt.xlabel('$\Delta \Delta G$ (kcal/mol)')
            plt.ylabel('number of reads/number of reads from input')
            plt.tight_layout()        
            """
            # generate a similar plot for random sequences
            fake_seqs = pd.Series(list(itertools.chain(*[[name]*val for name, val in motif_data.groupby('seq').size().iteritems()])),
                                  index=np.random.choice(random_data.index.tolist(), size=len(motif_data), replace=False))
    
    
            cols = ['input', 'rep1', 'rep2']
            random_data.loc[:, 'fake_seq'] = fake_seqs 
            counts = random_data.groupby('fake_seq')[cols].mean()*weights
            counts_bounds = pd.concat({(name, col):pd.Series(np.abs(bootstrap.ci(group.loc[:, col]*weights.loc[col], np.mean) - counts.loc[name, col]),
                                                             index=['eminus', 'eplus'])
                             for col in cols for name, group in random_data.groupby('fake_seq')})
            counts_err = pd.concat({name:group.unstack().loc[:, name].unstack()
                                    for name, group in counts_bounds.groupby(level=2)})
            mat = pd.concat([ddG, counts], axis=1)
            # plot
            plt.figure(figsize=(4,3))
            for col in cols:
                plt.scatter(mat.ddG, mat.loc[:, col])
                plt.errorbar(mat.ddG, mat.loc[:, col], yerr=[counts_err.loc[key].loc[mat.index.tolist(), col] for key in ['eminus', 'eplus']], fmt=',')
            plt.yscale('log')
            plt.xlabel('$\Delta \Delta G$ (kcal/mol)')
            plt.ylabel('norm. num. reads per site')
            plt.tight_layout()        
    
        if args.mode=='find_background_estimates':
            """Match each site with a site of similar expression and subtract."""
            logdata = np.log2(data.loc[:, ['input', 'tpm']]+0.05)
            fitparams = pd.concat({key:tectplots.fitGaussianKernel(logdata.loc[:, key], plot=True) for key in logdata}).unstack()
            
            zscores = (logdata - fitparams.center)/(fitparams.sigma)
            fillna = zscores.groupby('input')['tpm'].median()
            zscores.loc[zscores.tpm.isnull(), 'tpm'] = fillna.loc[zscores.loc[zscores.tpm.isnull()].input].values
            zscores.dropna(inplace=True)
            # find knn
            interp_info = clustering.cluster_knn_general(zscores.loc['random', ['tpm']], zscores.loc['motif', ['tpm']], k=100, threshold=0.5)
            
            # estimate counts from non specific reads per site.
            index_dict = {key: pd.MultiIndex.from_tuples([('random', idx) for idx in neighbors]) for key, neighbors in interp_info.neighbors.iteritems()}
            background_estimates = {}
            for i, (motif, neighbors) in enumerate(interp_info.neighbors.iteritems()):
                if (i%int(len(interp_info)/20.))==0:
                    print 'completed %4.1f%%'%(i/(len(interp_info)/100.))
                
                background_estimates[motif]  = np.exp(np.log(data.loc[index_dict[motif], cols]+1).mean())-1    
                
                #mat = data.loc['random'].loc[neighbors, cols].copy()
                #background_estimates[motif] = pd.concat({'geom_mean':np.exp(np.log(mat+1).mean())-1,'med':mat.median(), 'meen':mat.mean()}, names=['agg_method', 'sample'])
            background_estimates = pd.concat(background_estimates, names=['site'])
            expected_counts = background_estimates.unstack()
    
            count_sub = (motif_data.loc[:, cols] - expected_counts)
            
            norm_count_sub = (motif_data.loc[:, cols] - expected_counts)*weights
            norm_enrichment = ((norm_count_sub.rep1 + norm_count_sub.rep2)/(2*norm_count_sub.input))
            norm_enrichment.loc[motif_data.input==0] = np.nan
            
            norm_enrichment.replace([np.inf, -np.inf], np.nan, inplace=True)
            motif_data.loc[:, 'fold_enrichment_sub'] = norm_enrichment.astype(float)
            for col in ['fold_enrichment', 'fold_enrichment_sub']:
                mat = motif_data.groupby('seq')[['ddG', col]].median()
                plt.figure(figsize=(3,3))
                plt.scatter(mat.ddG, mat.loc[:, col])
                plt.title(col)
                plt.tight_layout()
                #plt.yscale('log')
            
        
        if args.mode=='compare_transcripts_with_binding_sites':
            """Find those genes with 'good' sites, and compare their GO terms +/- clip signal."""
            gene_clip_max = bed_data.loc[bed_data.binned_ddG==0].groupby('gene')['clip_enrichment'].max()
            genes_sorted = gene_clip_max.dropna().sort_values(ascending=False).index.tolist()
            
            go_terms = pd.read_table('GO_process.txt', index_col=0)
            gene_go_terms = {}
            for idx, row in go_terms.iterrows():
                
                gene_list = [a[1:a.find('  -  ')] for a in row.Genes.split(',') if a.find('  -  ')>=0]
                qvalue = pd.Series(-np.log10(row.loc['FDR q-value']), index=gene_list).rename('log10q')
                gene_go_terms[idx] = qvalue
            gene_go_terms = pd.concat(gene_go_terms, axis=1).fillna(0)
            
            # cluster
            z = sch.linkage(gene_go_terms, method='ward')
            z_row = sch.linkage(gene_go_terms.transpose(), method='ward')
    
            col_clusters = pd.Series(sch.fcluster(z_row, 6, criterion='maxclust'), index=gene_go_terms.columns)
            cluster_colors = pd.Series({i:c for c, i in zip(sns.color_palette('Set1', n_colors=6), np.unique(col_clusters))})
            col_colors = [cluster_colors.loc[i] for i in col_clusters]
            sns.clustermap( gene_go_terms, row_linkage=z, col_linkage=z_row, col_colors=col_colors, yticklabels=False, xticklabels=False)
            
            # reduce by col cluster
            gene_clustered = pd.concat([gene_go_terms.transpose(), col_clusters.rename('cluster')], axis=1).groupby('cluster').median().transpose()
            index_row = (gene_clustered>0).any(axis=1)
            z_row = sch.linkage(gene_clustered.loc[index_row], method='ward')
            sns.clustermap( gene_clustered.loc[index_row], row_linkage=z_row, col_cluster=False, col_colors=cluster_colors.tolist(), yticklabels=False, xticklabels=False, vmin=0, vmax=4)
            sys.exit()
            pd.concat([go_terms.loc[:, ['Description', 'FDR q-value']], col_clusters.rename('cluster')], axis=1).sort_values(['cluster', 'FDR q-value'])
            
            # plot reduced clusters
            index_col = col_clusters!=1
            index_row = (gene_go_terms.loc[:, index_col]>0).sum(axis=1)>10
            z_row = sch.linkage(gene_go_terms.loc[index_row], method='ward')
            z_col = sch.linkage(gene_go_terms.transpose().loc[index_col], method='ward')
            col_clusters = pd.Series(sch.fcluster(z_col, 5, criterion='maxclust'), index=gene_go_terms.loc[:, index_col].columns)
            col_colors = [cluster_colors.loc[i] for i in col_clusters]
            sns.clustermap( gene_go_terms.loc[index_row, index_col], row_linkage=z_row, col_linkage=z_col, col_colors=col_colors)
            #sns.clustermap( gene_go_terms.loc[index_row, index_col], )
        
        if args.mode=='compare_transcripts_more':
            """Find genes with good binding sites, but no enrichmnet, or lots of enrichmnent, and compare to all genes."""
            transcribed_genes = pd.read_table('RNAseq/transcript_quant/rna_seq_combined.tpm.above_0.01_both.dat').transcript_id
            biomart_converter = pd.read_table('annotations/ensemble_gene_converter_biomart.txt', header=0, names=['gene_id', 'transcript_id', 'gene_name', 'refseq_mRNA', 'refseq_ncRNA'])
            refseq_transcripts = biomart_converter.set_index('transcript_id').loc[[s.split('.')[0] for s in transcribed_genes if s.find('ENST')==0]].refseq_mRNA.dropna().unique()
            
            # find transcripts with no clip sites
            gene_clip_max = bed_data.loc[bed_data.binned_ddG==1].groupby('gene')['clip_enrichment'].max()
            binned_clip = pd.Series(np.digitize(gene_clip_max.dropna(), [-10, 1, 4, 11]), index=gene_clip_max.dropna().index)
            genes_no_clip_sites = [s for s, b in zip(refseq_transcripts, np.in1d(refseq_transcripts, gene_clip_max.index.tolist())) if not b]
            genes_low_enrichment = binned_clip.loc[binned_clip==1].index.tolist()
            genes_med_enrichment = binned_clip.loc[binned_clip==2].index.tolist()
            genes_high_enrichment = binned_clip.loc[binned_clip==3].index.tolist()
            
            pd.Series(genes_low_enrichment).to_csv('analysis/hPUM2_1muts/go_terms/genes.low_enrichment.consensus_site.txt', index=False)
            
            # compare these using gorilla Go term analysis, export as xls and save as txt
            mydir = 'analysis/hPUM2_1muts/go_terms/'
            filenames = {'consensus_sorted':mydir + 'GO_process_consensus_sorted.txt',
                         'high_to_low':     mydir + 'GO_process_high_to_low.txt',
                         'high_to_no_site': mydir + 'GO_process_high_to_no_site.txt',
                         'low_to_no_site':  mydir + 'GO_process_low_to_no_site.txt',}
            
            qvalues = {}
            for key, filename in filenames.items():
                mat = processing.load_go_file_with_errors(filename)
                qvalues[key] = -np.log10(mat.fdr_qvalue.astype(float))
            qvalues = pd.concat(qvalues, axis=1).fillna(0)
            
            grouped = (qvalues>0).groupby([col for col in qvalues])
            cluster = pd.concat([pd.Series('_'.join(['%d'%i for i in name]), index=group.index) for name, group in grouped]).rename('cluster')
            distances = pd.concat([clustering.get_distance_from_median_general(qvalues.loc[group.index]) for name, group in grouped]).rename('distance')      
            
            z = sch.linkage(qvalues , method='ward')
            sns.clustermap(qvalues, yticklabels=False, row_linkage=z)
            sns.clustermap(pd.DataFrame(col_clusters.loc[qvalues.index.tolist()]), yticklabels=False, row_linkage=z, cmap='Set1', col_cluster=False)   
                             
            
            
    if args.mode == 'compare_between_pos_5_variants':
        """Compare the consensus and similar sites"""
        bed_data = mut_bed_data.loc[mut_bed_data.same_strand].copy()
        bed_data.loc[:, 'seq'] =  pd.read_table('analysis/hPUM2_1muts/sequence/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.fasta.dat',
                                                header=None,  names=['name', 'seq'], index_col='name', squeeze=True)
        rc_seqs = []
        for name, group in bed_data.groupby('strand'):
            if name == '+':
                rc_seqs.append(group.seq.str.upper())
            else:
                rc_seqs.append(pd.Series({key: seqfun.rc(s.upper()) for key, s in group.seq.iteritems()}))
        rc_seqs = pd.concat(rc_seqs)
        bed_data.loc[:, 'seq_rc'] = rc_seqs
        
        # annotate seqs with position
        consensus = 'TGTATATA'
        index_list = {}
        for s in bed_data.seq_rc.unique():
            vec = [(i, s2) for i, (s1, s2) in enumerate(zip(consensus, s)) if s1!=s2]
            if len(vec)==1:
                index_list[s] = (vec[0])
            elif len(vec) == 0:
                index_list[s] = ((-1, ''))
            else:
                print 'double mutant found!'
                break
        index_list = pd.Series(index_list)
    
        # find the variants at position 5 (1 index)
        position = pd.Series({idx:vec[0] for idx, vec in index_list.iteritems()})
        base_sub = pd.Series({idx:vec[1] for idx, vec in index_list.iteritems()})
        bed_data.loc[:, 'position'] = position.loc[bed_data.seq_rc].values
        bed_data.loc[:, 'base_sub'] = base_sub.loc[bed_data.seq_rc].values
        
        # plot
        sns.factorplot(data=bed_data.loc[(bed_data.position==4)|(bed_data.position==-1)],
                       x='base_sub', y='clip_enrichment', kind='box', notch=True, palette=['0.5'], aspect=0.7); plt.ylim(-5, 10)
    
      
    if args.mode == 'compare_to_transcript_abundance':
        """Compare the clip enrichment to the transcript abundance."""
        transcript_quant_all = pd.read_table('RNAseq/transcript_quant/rna_seq_combined.tpm.above_0.01_both.dat')
        transcript_quant = transcript_quant_all.set_index('transcript_id').mean(axis=1)
        peak_id = pd.read_table('analysis/hPUM2_1muts/beds/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.transcript.peak_id',
                                header=None, names=['name', 'exon_id', 'transcripts'])
        
        # associate each motif with a single transcript based on the most expressed transcript that overlaps with the motif site
        motif_id = pd.Series({motif: transcript_quant.loc[transcript_str.split(',')].idxmax() for i, motif, transcript_str in peak_id.loc[:, ['name', 'transcripts']].itertuples()})
    
        bed_data = mut_bed_data.loc[mut_bed_data.same_strand].copy()
        bed_data.loc[:, 'tpm'] = transcript_quant.loc[motif_id.loc[bed_data.index.tolist()]].values
        bed_data.loc[:, 'log_tpm'] = np.log2(bed_data.tpm)
        
        # plot
        xlim = np.array([-5, 10])
        x = bed_data.loc[bed_data.score>11].dropna(subset=['log_tpm', 'clip_enrichment']).log_tpm
        y = bed_data.loc[bed_data.score>11].dropna(subset=['log_tpm', 'clip_enrichment']).clip_enrichment
        plt.figure(figsize=(3,3))
        plt.scatter(x, y, marker='.', c='0.5', edgecolor='b', linewidth=0.3)
        slope, intercept, rvalue, pvalue, stderr = st.linregress(x, y)
        plt.plot(xlim, xlim*slope + intercept, 'k', linewidth=1)
        plt.xlim(xlim); plt.ylim(xlim)
        sns.despine()
        plt.annotate('rsq=%s')
        
    if args.mode=='find_NORAD_signal':
        """Find sites within the NORAD ncRAN."""
        
        # load the unfiltered bed file.
        a = pd.read_table('analysis/hPUM2_1muts/beds/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.bed')
    
        
        num_norad = a.loc[a.loc[:, 'Gene Type']=='ncRNA'].groupby('Nearest Ensembl').size().loc['ENSG00000260032']
        data = a.loc[a.loc[:, 'Gene Type']=='ncRNA'].groupby(['Nearest Ensembl', 'Peak Score']).size().reset_index()
        order = [4.161568, 11.066319]
        g = sns.FacetGrid(data=data, col='Peak Score', sharey=False, sharex=False, col_order=order, hue='Peak Score', hue_order=order); g.map(sns.distplot, 0, kde=False)
        for i, ax in enumerate(g.axes.flat):
            ax.axvline(num_norad.loc[num_norad.loc[:, 'Peak Score']==order[i]].shape[0])
        
    if args.mode=='compare_nucleotide_content':
        """Given the nucleotide content around consensus sites, copmare to clip enrichmnet."""
        
        nuc_content = pd.read_table('analysis/hPUM2_0muts/sequence/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.nuc_content', index_col=0)
        bed_data = mut_bed_data.copy()
        
        cols = [ '9_num_A', '10_num_C', '11_num_G','12_num_T',]
        cols = [ '7_pct_at',
     '8_pct_gc',]
        data = pd.melt(pd.concat([bed_data.loc[bed_data.score>11].clip_enrichment, nuc_content.loc[:, cols]], axis=1), id_vars=['clip_enrichment'])
        g = sns.FacetGrid(data, col='variable'); g.map( sns.regplot, 'value', 'clip_enrichment', marker='.')
        
        
    if args.mode=='parse_cds_ends':
        """Save only cds ends of genes that have puf sites."""
        bed_file = pd.read_table('annotations/refseq/hg38_refGene.cds_end.st.bed', header=None, names=bed_fields)
        index = np.in1d(bed_file.name,  mut_bed_data.gene.unique())
        bed_file.loc[index].to_csv('analysis/hPUM2_1muts/beds/hg38_refGene.cds_end.st.subset_genes.bed', sep='\t', index=False, header=False)
    if args.mode=='parse_tx_ends':
        """Save only cds ends of genes that have puf sites."""
        bed_file = pd.read_table('annotations/refseq/hg38_refGene.tx_end.st.bed', header=None, names=bed_fields)
        index = np.in1d(bed_file.name,  mut_bed_data.gene.unique())
        bed_file.loc[index].to_csv('analysis/hPUM2_1muts/beds/hg38_refGene.tx_end.st.subset_genes.bed', sep='\t', index=False, header=False)
        
        
    if args.mode=='plot_distance_to_cds_ends':
        """load output of bedtools closest command."""
        # bedtools closest -d -a $bed -b analysis/hPUM2_1muts/beds/hg38_refGene.cds_end.st.subset_genes.bed | awk -F "\t" '{OFS="\t"}{if ($8==$13) print $4, $6, $13, $15, $16}' > analysis/hPUM2_1muts/location/distance_to_closest_cds_end.txt
        distance_mat = pd.read_table('analysis/hPUM2_1muts/location/distance_to_closest_tx_end.txt', header=None, names=['name', 'motif_strand', 'gene_name', 'gene_strand', 'distance'], index_col='name')
        data = pd.concat([mut_bed_data, np.log10(distance_mat.distance+1)], axis=1)
        g = sns.FacetGrid(data.loc[(data.score>11)&(data.same_strand)], ); g.map(plt.scatter, 'distance', 'clip_enrichment', marker='.')
        x, y = [data.loc[(data.score>11)&(data.same_strand)].dropna().loc[:, col] for col in ['distance', 'clip_enrichment']]
        slope, intercept, rvalue, pvalue, stderr = st.linregress(x, y)
        xlim = np.array([0, 6])
        plt.plot(xlim, xlim*slope+intercept, 'k--')
    
    if args.mode=='bin_by_number_of_reads':
        """Load clip output, find number of reads, and bin clip enrichment."""
        filename = 'analysis/hPUM2_1muts/clip_signal/hPUM2.1_muts.exons.st.merge_transcript.above_0.01_both.st.ann.filt.clip.counts'
        data_agg2 = pd.read_csv(filename, index_col=0)
        num_reads = (data_agg2/data_agg2.replace(0, np.nan).min())
        
        # from --mode compare_between_pos_5_variants
        num_reads_signal_log = np.log2((num_reads.rep1 + num_reads.rep2 + 1)/2.)
        num_reads_input_log = np.log2(num_reads.input + 1)
    
        index_consensus_plus = bed_data.loc[((bed_data.position==4)|(bed_data.position==-1))&(bed_data.base_sub!='G')].dropna(subset=['clip_enrichment']).index.tolist()
        
        sorted_signal =  num_reads_input_log.loc[index_consensus_plus].replace(0, np.nan).dropna().sort_values()
        num_in_bin = len(index_consensus_plus)/7
        index_low = sorted_signal[num_in_bin:num_in_bin*2].index.tolist()
        index_high = sorted_signal[-num_in_bin*2:-num_in_bin].index.tolist()
        plt.figure(figsize=(4,3))
        for i, index in enumerate([index_low, index_high]):
            tectplots.distplot_kde_scaled(bed_data.loc[index].clip_enrichment, color=sns.color_palette()[i])
        num_bins = 7
        binned_signal = pd.Series(np.digitize(sorted_signal, bins=sorted_signal.quantile(np.linspace(0, 1, num_bins+1))), index=sorted_signal.index)
        
        # combine and plot
        a = pd.concat([num_reads_signal_log.rename('signal_counts_log'),
                       num_reads_input_log.rename('input_counts_log'),
                       bed_data.clip_enrichment, binned_signal.rename('binned_counts')], axis=1).loc[index_consensus_plus]
        g = sns.FacetGrid(data=pd.melt(a, id_vars=['clip_enrichment', 'binned_counts']), col='variable'); g.map(plt.scatter, 'value', 'clip_enrichment', marker='.', s=4)
        
        # bin and plot
    
        num_measurements = a.groupby('binned_counts').size()
        std_dev = a.groupby('binned_counts')['clip_enrichment'].std().rename('clip_enrichment_std')
        #std_dev_bounds = pd.concat([pd.Series(bootstrap.ci(group.clip_enrichment, np.std, method='pi'), index=['lb', 'ub']).rename(name)
        #                            for name, group in a.groupby('binned_counts')], axis=1).transpose()
        #std_dev_err = (std_dev_bounds.ub - std_dev_bounds.lb)/2.
        std_dev_err = std_dev/np.sqrt(num_measurements)
        x_values = pd.concat([binned_signal.rename('bin_name'), sorted_signal.rename('signal')], axis=1).groupby('bin_name')['signal'].median()
        (pd.concat([std_dev, x_values.rename('input_signal'), std_dev_err.rename('std_err') ], axis=1).
         plot(x='input_signal', y='clip_enrichment_std', yerr='std_err', marker='o', xlim=[5, 11], figsize=(3,3)))
        plt.xlabel('log2 (input # reads)')
        plt.ylabel('std. dev. of log2 clip enrichment')
        plt.tight_layout()
    
    if args.mode=='test_pwms':
        """Load the pwm for single mutants, and determine the right score"""
        from hjh import mutations
        
        consensus_seq = 'UGUAUAUAG'
        Single = predictions.Single()
        eps = 0.1
        max_ddG = 4 # kcal/mol
        # make up to 3 mutations in consensus_seq
        seqs = [consensus_seq]
        for i in range(3):
            seqs = set(list(itertools.chain(*[mutations.singles(s) for s in seqs])))
        seqs = list(seqs)   
        
        ddGs = pd.Series([Single.get_ddG(s) for s in seqs], index=seqs).rename('ddG')
        pwm2 = Single.pwm.copy()
        pwm2.loc[8] = [0.25]*4
        scores2 =  pd.Series([predictions.log_odds_score(s, Single.pwm) for s in seqs], index=seqs).rename('scores_old')
        scores =  pd.Series([predictions.log_odds_score(s, pwm2) for s in seqs], index=seqs).rename('scores')

        data = pd.concat([ddGs, scores, scores2], axis=1).reset_index().rename(columns={'index':'seq'})
        data = pd.concat([data, pd.concat({i:vec for i, vec in enumerate(list(data.seq.str))}, axis=1)], axis=1)
        data.loc[:, 'above_threshold'] = data.scores > 2
        data.loc[:, 'above_threshold_old'] = data.scores_old > 2
        g = sns.FacetGrid(data, hue=7); g.map(plt.scatter, 'ddG', 'scores', marker='.')
        scores_at_threshold = data.loc[np.abs(data.ddG - max_ddG) < eps].scores
        score_threshold = scores_at_threshold.median() - 2*scores_at_threshold.std()
        plt.axvline(max_ddG, linewidth=1, linestyle=':'); plt.axhline(score_threshold, linewidth=1, linestyle='--')
        
        sys.exit()
        # for those above the score threhsold, evaluate how many duokicates in the set of insertions
        seq_subset =  data.loc[(data.ddG < max_ddG)].seq.tolist()
        
        insertion_penalties_ddG = pd.read_table('annotations/RNAmap/insertion_penalties.dat', index_col=0)
        insertion_penalties_kd = np.exp(insertion_penalties_ddG/predictions.RT).replace(np.nan, np.inf)

        # find PWMs for single insertions
        pwm_ins3 = pd.concat({i:predictions.insert_base_in_pwm(i, Single.pwm, rel_kds=insertion_penalties_kd.loc[i]) for i in [3]}).groupby(level=1).mean()
        pwm_ins45 = pd.concat({i:predictions.insert_base_in_pwm(i, Single.pwm, rel_kds=insertion_penalties_kd.loc[i]) for i in [4,5]}).groupby(level=1).mean()
        pwm_ins6 = pd.concat({i:predictions.insert_base_in_pwm(i, Single.pwm, rel_kds=insertion_penalties_kd.loc[i]) for i in [6]}).groupby(level=1).mean()

        seqs_inserted = {}
        ddGs_ins = {}
        scores_ins = {}
        for position, penalties in insertion_penalties_ddG.dropna(how='all').iterrows():
            for s in seq_subset:
                for base, penalty in penalties.dropna().iteritems():
                    seq = s[:position] + base + s[position:]
                    key = (s, position, base)
                    seqs_inserted[key] = seq
                    ddGs_ins[key] = ddGs.loc[s] + penalty # its the single base sub penalty, + penalty, assuming there isn't a better register
                    if position == 3:
                        scores_ins[key] = predictions.log_odds_score(seq, pwm_ins3)
                    elif position == 4 or position == 5:
                        scores_ins[key] = predictions.log_odds_score(seq, pwm_ins45)
                    elif position == 6:
                        scores_ins[key] = predictions.log_odds_score(seq, pwm_ins6)
        ddGs_ins = pd.Series(ddGs_ins).rename('ddG')
        seqs_inserted = pd.Series(seqs_inserted).rename('seq')
        scores_ins = pd.Series(scores_ins).rename('score')
        data_ins = pd.concat([ddGs_ins, scores_ins, seqs_inserted], axis=1)
        for name, group in data_ins.groupby(level=1):
            scores_at_threshold = group.loc[np.abs(group.ddG - max_ddG) < eps].score
            score_threshold = scores_at_threshold.median() - 2*scores_at_threshold.std()
            g = sns.FacetGrid(group); g.map(plt.scatter, 'ddG', 'score', marker='.')
            plt.axvline(max_ddG, linewidth=1, linestyle=':'); plt.axhline(score_threshold, linewidth=1, linestyle='--')
            # decided to just keep score threshold as it is for single base subs, as it is very accomodating
            
        # double insertions
        double_insertion_penalties_ddG = pd.read_table('annotations/RNAmap/double_insertion_penalties.dat', index_col=0)
        pwm_ins2_3 = pd.concat({i:predictions.insert_base_in_pwm(i, Single.pwm, num_insertions=2) for i in [3]}).groupby(level=1).mean()
        pwm_ins2_45 = pd.concat({i:predictions.insert_base_in_pwm(i, Single.pwm, num_insertions=2) for i in [4,5]}).groupby(level=1).mean()

        # save pwms
        outdir = 'annotations/RNAmap/motifs/'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        length = pwm_ins2_3.shape[0]
        #odds_threshold = 1.62
        for name, pwm, odds_threshold in zip(['hPUM2', 'hPUM2_flip3', 'hPUM2_flip45', 'hPUM2_flip6', 'hPUM2_2flip3', 'hPUM2_2flip45'],
                                             [Single.pwm, pwm_ins3, pwm_ins45, pwm_ins6, pwm_ins2_3, pwm_ins2_45],
                                             [2, 4, 4, 4, 5, 5]):
            pwm = pwm.copy()
            # make all the same length by appending to the 3' end
            pwm = predictions.insert_base_in_pwm(len(pwm), pwm, num_insertions=length - len(pwm))
            # save
            filename = os.path.join(outdir, '%s.motif'%name)
            predictions.save_pwm(pwm, filename, odds_threshold, name=name)
        