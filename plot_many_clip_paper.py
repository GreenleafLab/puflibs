#!/usr/bin/env python
"""
This script may be used to plot all the plots in the paper figures.

Run:
python plot_many_clip_paper.py --mode plot_clip_footprint --data $datafilename 
python plot_many_clip_paper.py --mode plot_clip_vs_ddG --data $datafilename 
python plot_many_clip_paper.py --mode plot_overestimation_with_noflip_model --data $datafilename 
python plot_many_clip_paper.py --mode plot_annotation_vs_ddG --data $datafilename 
python plot_many_clip_paper.py --mode compare_consensus_sites_for_ss_structure --data $datafilename 
"""

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
from puflibs import processing, predictions, seqmodel, variables



# import args
parser = argparse.ArgumentParser()
parser.add_argument('--mode', help='which analysis to run')
parser.add_argument('--data', help='name of final data file. default = "analysis/output/analysis/output/all_unprocessed_st_merged.00.hPUM2_all.random_5e+06.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.temp37.combined_data.01.02.03.04.05.06.07.08.09.combined_data.gz',
                    default='analysis/output/all_unprocessed_st_merged.00.hPUM2_all.random_5e+06.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.temp37.combined_data.01.02.03.04.05.06.07.08.09.combined_data.gz')
parser.add_argument('--ss_data', help='name of file with the secondary structure ddG for'
                    'consensus sites. default = "analysis/output/hPUM2_all.first_register_consensus.temp37.dG.combined_data.gz',
                    default='analysis/output/hPUM2_all.first_register_consensus.temp37.dG.combined_data.gz')

parser.add_argument('-t', '--temperature', type=float, help='temperature in celcius. default = 37',
                    default=37.)
    
if __name__ == '__main__':

    args = parser.parse_args()
    temperature = args.temperature    
    data = pd.read_table(args.data, compression='gzip')
    


    
    if args.mode == 'compare_rna_seq':
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
        
        # get NORAD
        tpm_data.loc['ENST00000565493']

    elif args.mode == 'find_consensus_sites_for_ss_structure':
        """To evaluate secondary structure effects, we want consensus sites in one register"""        
        filenames = subprocess.check_output('find analysis/effects/temp_%d -mindepth 2 -name "*affinity.gz" | sort'%temperature, shell=True).strip().split()
        kT = seqmodel.get_ddG_conversion(temperature)
        data_subset = {}
        for i, filename in enumerate(filenames):
            split_id = int(filename.split('/')[3].split('_')[-1])
            if split_id in range(10): # exclude split_10 which is peaks
                seqeffect = pd.read_table(filename, compression='gzip', index_col=0)
                name_list = seqeffect.loc[seqeffect.noflip_0 < 0.5].index.tolist()
                
                data_subset[split_id] = (data.groupby('split_id').get_group(split_id).drop('split_id', axis=1).set_index('name').
                                         loc[name_list].copy())
            
        data_subset = pd.concat(data_subset, names=['split_id', 'name2']).dropna(subset=['chrm', 'start', 'stop']).reset_index()
        data_subset.loc[:, 'name'] = data_subset.split_id.astype(str) + ',' + data_subset.name2
        for col in ['start', 'stop']:
            data_subset.loc[:, col] = data_subset.loc[:, col].astype(int)
         
        #data_subset.loc[:, variables.bed_fields].to_csv('analysis/beds/hPUM2_all.first_register_consensus.bed',
        #                                                sep='\t', index=False, header=False, float_format='%.2f')
    
    elif args.mode == 'compare_consensus_sites_for_ss_structure':
        """Evaluate secondary structure effects for 8nt constraint and 11nt constraint"""
        filename = args.ss_data # analysis/output/hPUM2_all.first_register_consensus.temp37.dG.combined_data.gz
        ss_ddG_table = pd.read_table(filename, compression='gzip', index_col=(0,1))                       


        # find clip signal
        data_subset = data.set_index(['split_id', 'name']).loc[ss_ddG_table.index.tolist()].copy()
        med_background_val = data.loc[(data.tpm>0)&(data.ddG > 4.5)].clip_signal_per_tpm.median()
        clip_signal_per_tpm_fold = data_subset.clip_signal_per_tpm/med_background_val
        
        # bin the ss_ddGs
        binedges = ss_ddG_table.stack().quantile(np.linspace(0, 1, 15))
        ss_ddG_binned = pd.cut(ss_ddG_table.stack(), bins=binedges, precision=1, include_lowest=True).rename('ss_ddG')
        order = pd.DataFrame(ss_ddG_binned).groupby('ss_ddG').first().index.tolist()
        
        # find entries in the first bin in at least one energy category
        first_bin = ss_ddG_binned.loc[ss_ddG_binned == order[0]].unstack().index.tolist()
        ymax = clip_signal_per_tpm_fold.loc[first_bin].median()

        # plotting function
        func = functools.partial(sns.factorplot, x='ss_ddG', y='clip_signal_per_tpm',  estimator=np.median,
                                         errwidth=0.5,capsize=1, linestyles='', marker='.')
        ylim = [0.5, 500]       

        # plot with equally spaced bins
        binedges = np.hstack([np.arange(0, 5, 0.5), np.logspace(np.log10(4.5), np.log10(ss_ddG_table.stack().max()), 3)])
        ss_ddG_binned = pd.cut(ss_ddG_table.stack(), bins=binedges, precision=1, include_lowest=True).rename('ss_ddG')
        order = pd.DataFrame(ss_ddG_binned).groupby('ss_ddG').first().index.tolist()
        x_order = (binedges[1:] + binedges[:-1])*0.5

        
        count_cutoff = 20
        all_data = {}
        for i, window_size in enumerate([10, 20, 40, 80, 160]):
            data_toplot = (pd.concat([clip_signal_per_tpm_fold,
                                      ss_ddG_binned.unstack().loc[:, 'ss_ddG_%d'%window_size].rename('ss_ddG')], axis=1)
                           .loc[(data_subset.tpm>0)&(data_subset.in_transcript)])
            to_plot_categories = [category for category, count in data_toplot.ss_ddG.value_counts().iteritems() if count > count_cutoff]
            
            all_data[window_size] = data_toplot.loc[np.in1d(data_toplot.ss_ddG, to_plot_categories)]
            #all_data[window_size] = data_toplot
        all_data = pd.concat(all_data, names=['window_size', 'split_id', 'name']).reset_index()

        g = func(data=all_data, hue='window_size', hue_order=[20, 40, 160])
        plt.yscale('log')
        plt.xticks(rotation=90)
        plt.axhline(1, color='0.5', linestyle=':')
        
        # plot expectation line
        y = ymax*np.exp((x_order-x_order[0])/seqmodel.get_ddG_conversion(temperature))
        plt.plot(np.arange(len(order)), y, 'k--')
        plt.ylim(ylim)
            
        
            


    
    elif args.mode == 'plot_clip_footprint':
        """Load the clip signal counts and plot the aggregate footprint around consensus sites"""
        
        filenames_rep1 = subprocess.check_output('find analysis/clip -mindepth 2 -maxdepth 2 -name "*tracks.txt.gz" | grep split | grep hPUM2_all | grep -v input | grep -v rep2 | sort', shell=True).strip().split()
        filenames_rep2 = subprocess.check_output('find analysis/clip -mindepth 2 -maxdepth 2 -name "*tracks.txt.gz" | grep split | grep hPUM2_all | grep -v input | grep -v rep1 | sort', shell=True).strip().split()
        filenames_input = subprocess.check_output('find analysis/clip -mindepth 2 -maxdepth 2 -name "*tracks.txt.gz" | grep split | grep hPUM2_all | grep -v rep2 | grep -v rep1 | sort', shell=True).strip().split()

        # load counts
        
        footprints = {}
        for i, filename in enumerate(filenames_rep1):
            split_id = int(filename.split('/')[2].split('_')[-1])
            print split_id
            footprint = pd.read_csv(filename, compression='gzip', index_col=0)
            subset_consensus = data.loc[
                (data.split_id == split_id)&
                (data.ddG < 0.5)&
                (data.tpm > 0)&
                (data.in_transcript)].name
            footprints[split_id] = footprint.loc[subset_consensus]
            
        footprints = pd.concat(footprints, names=['split_id', 'name'])
        
        # plot
        interval_radius=40
        offset=15
        xvalues = np.arange(-250, 251)
        
        plt.figure(figsize=(3,3))
        plt.plot(xvalues, footprints.mean())
        plt.axvline(-interval_radius-offset, color='k', linestyle=':')
        plt.axvline(+interval_radius-offset, color='k', linestyle=':')
        plt.ylim(0, 0.6)
        plt.xlim(-225, 225)
        plt.savefig('lineplot.footprint.pdf')
    

    
    elif (args.mode == 'plot_flip_annot' or args.mode == 'plot_clip_vs_ddG' or
          args.mode == 'plot_annotation_vs_ddG' or
          args.mode == 'plot_overestimation_with_noflip_model' or
          args.mode == 'plot_intron_clip_subset' or
          args.mode=='plot_ss_ddG_versus_clip' or args.mode =='plot_versus_proximity'):
        """Using the New model to find effects, find ddG"""
        pass
        # annotatie flipped/not flipped
        data.loc[:, 'flip_annot'] = 'noflip'
        data.loc[data.ddG_noflip - data.ddG > 0.5, 'flip_annot'] = 'flip'
        
        # find fold enrichment above expected bacground
        med_background_val = data.loc[(data.tpm>0)&(data.ddG > 4.5)].clip_signal_per_tpm.median()
        med_background_val_input = data.loc[(data.tpm>0)&(data.ddG > 4.5)].clip_input_per_tpm.median()
        data.loc[:, 'clip_signal_per_tpm_fold'] = data.clip_signal_per_tpm/med_background_val
        data.loc[:, 'clip_input_per_tpm_fold'] = data.clip_input_per_tpm/med_background_val_input

        # determine the subset of data with expression and within a transcript
        subdata = data.loc[(data.tpm > 0)&data.in_transcript]

        
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
            bins = np.linspace(0, 5)
            plt.figure(figsize=(3,3)); sns.distplot(data.loc[(data.tpm>0)&(~data.is_random)].dddG_noflip, color=sns.color_palette()[0], bins=bins, kde=False)
                

        elif args.mode == 'plot_clip_vs_ddG' or args.mode=='plot_versus_proximity':

            """bin by ddG and plot"""
            ddG_binedges = np.hstack([np.linspace(data.ddG.min(), 4.5, 25), data.ddG.max()])
            subdata.loc[:, 'binned_ddG'] = pd.cut(subdata.ddG, ddG_binedges,
                                                  include_lowest=True, precision=2)
            subdata.loc[:, 'binned_logtpm'] = pd.cut(np.log10(subdata.tpm), (np.log10(subdata.tpm).replace(-np.inf, np.nan).dropna()).quantile(np.linspace(0, 1, 4)), include_lowest=True)
            subdata.loc[:, 'binned_min_dist'] = pd.cut(subdata.min_bases_to_TGTA, [0, 50, 100, 300], include_lowest=True, precision=0)
            order = subdata.groupby('binned_ddG').first().index.tolist()
            
            # plot
            ylim = [0.5, 100]
            ymax = subdata.groupby('binned_ddG').get_group(order[0])['clip_signal_per_tpm_fold'].median() 
            x_order = pd.Series({name:0.5*(group.ddG.max() + group.ddG.min()) for name, group in subdata.groupby('binned_ddG')})
            
            func = functools.partial(sns.factorplot, data=subdata, estimator=np.median,
                                         errwidth=0.5,capsize=1, linestyles='', marker='.')
            
            # plot signal, colored by flip/noflip
            for yval in ['clip_input_per_tpm_fold', 'clip_signal_per_tpm_fold']:
            
                g = func(x='binned_ddG', y=yval, hue='flip_annot', hue_order=['noflip', 'flip']);
                plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
                plt.axhline(1, color='0.5', linestyle='--');
                plt.yscale('log')
                plt.ylim(ylim)
                # plot expected line
                y = ymax*np.exp(x_order/seqmodel.get_ddG_conversion(temperature))
                plt.plot(np.arange(len(order)), y.loc[order], 'k--')
                plt.savefig('scatterplot.%s.vs.binned_ddG.pdf'%yval)
            
            # plot the numbers of flipped sites per bin
            subdata.groupby(['binned_ddG', 'flip_annot']).size().unstack().loc[:, ['noflip', 'flip']].plot(kind='bar', stacked=True, width=0.8, figsize=(3,3));
            plt.savefig('barplot.num_flipped.binned_ddG.pdf')


            # plot signal, colored by annotation
            annotation_order = ["3' UTR", "exon"]
            annotation_colors = ['#f7931d', '0.7']
            for yval in ['clip_input_per_tpm_fold', 'clip_signal_per_tpm_fold']:
            
                g = func(x='binned_ddG', y=yval, hue='annotation', hue_order=annotation_order, palette=annotation_colors);
                plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
                plt.axhline(1, color='0.5', linestyle='--');
                plt.yscale('log')
                plt.ylim(ylim)
                # plot expected line
                y = ymax*np.exp(x_order/seqmodel.get_ddG_conversion(temperature))
                plt.plot(np.arange(len(order)), y.loc[order], 'k--')
                plt.savefig('scatterplot.%s.vs.binned_ddG.by_annotation.pdf'%yval)

            # plot by annotation again but only plot sites at least 100 nt away
            annotation_order = ["3' UTR", "exon"]
            annotation_colors = ['#f7931d', '0.7']
            func2 = functools.partial(sns.factorplot, data=subdata.loc[subdata.binned_min_dist=='(100, 300]'], estimator=np.median,
                                         errwidth=0.5,capsize=1, linestyles='', marker='.')
            for yval in ['clip_input_per_tpm_fold', 'clip_signal_per_tpm_fold']:
            
                g = func2(x='binned_ddG', y=yval, hue='annotation', hue_order=annotation_order, palette=annotation_colors);
                plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
                plt.axhline(1, color='0.5', linestyle='--');
                plt.yscale('log')
                plt.ylim(ylim)
                # plot expected line
                y = ymax*np.exp(x_order/seqmodel.get_ddG_conversion(temperature))
                plt.plot(np.arange(len(order)), y.loc[order], 'k--')
                plt.savefig('scatterplot.%s.vs.binned_ddG.by_annotation.greater_than_100nt.pdf'%yval)
            
    
            # plot signal, colored by distance to nearest site
            min_dist_order = ['[0, 50]', '(50, 100]', '(100, 300]']
            min_dist_colors = ['#f6935a', '#ab665c', '#603e4a']
            for yval in ['clip_input_per_tpm_fold', 'clip_signal_per_tpm_fold']:
            
                g = func(x='binned_ddG', y=yval, hue='binned_min_dist', hue_order=min_dist_order,
                         palette=min_dist_colors)
                plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
                plt.axhline(1, color='0.5', linestyle='--');
                plt.yscale('log')
                plt.ylim(ylim)
                # plot expected line
                y = ymax*np.exp(x_order/seqmodel.get_ddG_conversion(temperature))
                plt.plot(np.arange(len(order)), y.loc[order], 'k--')
                plt.savefig('scatterplot.%s.vs.binned_ddG.by_mindist.pdf'%yval)

            # plot fraction of annotations in distance bins
            subdata_stable  = subdata.loc[(subdata.ddG < 2)]
            num_in_annotations = subdata_stable.groupby(['annotation', 'binned_min_dist']).size().unstack()
            (num_in_annotations.transpose()/num_in_annotations.sum(axis=1)).loc[min_dist_order,annotation_order].transpose().plot(kind='bar', stacked=True, colors=min_dist_colors)
            plt.savefig('barplot.num_annotations.binned_min_dist.pdf')
            

            

            
        # plot by type
        elif args.mode == 'plot_rna_expression_vs_ddG':
            """Compare to RNA expression."""
            subdata = data.loc[(data.tpm>0)&(data.binned_ddG==order[0])].copy()
            subdata.loc[:, 'binned_tpm'] = np.digitize(subdata.tpm, np.logspace(-2, 3, 10))
            # plot the effect of expression
            # plot signal, colored by flip/noflip
            tpm_order = subdata.groupby('binned_tpm').first().index.tolist()
            for yval in ['clip_input_per_tpm_fold', 'clip_signal_per_tpm_fold']:
            
                g = func2(x='binned_ddG', y=yval, hue='binned_tpm', );
                plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
                plt.axhline(1, color='0.5', linestyle='--');
                plt.yscale('log')
                plt.ylim(ylim)
                # plot expected line
                y = ymax*np.exp(x_order/seqmodel.get_ddG_conversion(temperature))
                plt.plot(np.arange(len(order)), y.loc[order], 'k--')
                plt.savefig('scatterplot.%s.vs.binned_ddG.pdf'%yval)            
            
            
        
        elif args.mode == 'plot_overestimation_with_noflip_model':
            """Try to see if clip signal is overestimated when using ddG_no flip model"""
            subdata = data.loc[(data.tpm>0)&data.in_transcript&(subdata.min_bases_to_TGTA>=100)]
            
            binedges = np.hstack([np.arange(-0.5, 5, 0.5), subdata.ddG.max()])
            binedges = np.hstack([subdata.loc[(subdata.ddG < 4.5)&(subdata.ddG_noflip < 4.5), ['ddG', 'ddG_noflip']].stack().quantile(np.linspace(0, 1, 20)).values,
                                  max(subdata.ddG.max(), subdata.ddG_noflip.max())])
            binedges = np.hstack([np.linspace(data.ddG.min(), 4.5, 25), data.ddG.max()])

            subdata.loc[:, 'binned_ddG'] = pd.cut(subdata.ddG, binedges,
                                                  precision=1, include_lowest=True)
            subdata.loc[:, 'binned_ddGnoflip'] = pd.cut(subdata.ddG_noflip, binedges,
                                                  precision=1, include_lowest=True)
            order = subdata.groupby('binned_ddG').first().index.tolist()
            x_order = (binedges[1:] + binedges[:-1])*0.5
            ymax = subdata.groupby('binned_ddG').get_group(order[0]).clip_signal_per_tpm_fold.median()
            ylim = [0.1, 100]
            func = functools.partial(sns.factorplot, estimator=np.median,
                                         errwidth=0.5,capsize=0.5, linestyles='', marker='.')
            
             
            
            #plot
            yval = 'clip_signal_per_tpm_fold'
            count_cutoff = 25
            for xval in ['binned_ddG', 'binned_ddGnoflip']:
                subsubdata = pd.concat([group for name, group in subdata.groupby([xval, 'flip_annot']) if len(group) > count_cutoff])
                g = func(data=pd.concat([subsubdata.loc[:,[xval, 'flip_annot']], subsubdata[yval] + 0.1], axis=1), x=xval, y=yval, hue='flip_annot', hue_order=['noflip', 'flip'], order=order);
                
                plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
                plt.axhline(1, color='0.5', linestyle='--');
                plt.yscale('log')
                plt.ylim(ylim)
                # plot expected line
                y = ymax*np.exp((x_order-x_order[0])/seqmodel.get_ddG_conversion(temperature))
                plt.plot(np.arange(len(order)), y, 'k--')
                plt.savefig('scatterplot.clip_signal_vs_%s.big_bins.greaterthan_100nt.pdf'%xval)
        
        elif args.mode == 'plot_intron_clip_subset':
            """Using a small subset of sites for which the clip signal was determined for intronic sites, find the enrichment"""
            data_intron = pd.read_table('analysis/output/split_0/all_unprocessed_st_merged.09.hPUM2_odds9.intron_filt.input.ENCFF786ZZB.R2.500.rep2.ENCFF732EQX.rep1.ENCFF231WHF.temp0.combined_data.gz')
            data_intron.loc[:, 'min_bases_to_TGTA'] = data_intron.loc[:, ['upstream_bases_to_TGTA', 'downstream_bases_to_TGTA']].min(axis=1)

            subdata_intron = data_intron.loc[(data_intron.tpm>0)&(data_intron.min_bases_to_TGTA>=100)]
            subdata_intron.loc[:, 'clip_signal_per_tpm_fold'] = subdata_intron.clip_signal_per_tpm/med_background_val
            subdata = data.loc[data.in_transcript&(data.tpm>0)&(data.min_bases_to_TGTA>=100)]
            
            # plot
            annotation_order = ["3' UTR", "exon"]
            annotation_colors = ['#f7931d', '0.7']
            g = sns.FacetGrid(data=subdata.loc[subdata.ddG<0.5], hue='annotation',
                              hue_order=annotation_order, palette=annotation_colors);
            g.map(tectplots.plot_cdf, 'clip_signal_per_tpm_fold');
            tectplots.plot_cdf(subdata_intron['clip_signal_per_tpm_fold'], color='#494979'); plt.xscale('log')

            

            
        elif args.mode == 'plot_annotation_vs_ddG':
            """Using New model to find effects, se how enrichment for UTR changes"""

            annotation_order = ["3' UTR", "exon", "5' UTR"]
            annotation_colors = ['#f7931d', '0.7', '#be1e2d']
            
            # examine only sites with some hPUM annotation
            subdata = data.loc[data.in_transcript&(data.name.str.find('hPUM') == 0)]
            
            # find fraction per ddG bin
            subdata.loc[:, 'binned_ddG'] = pd.cut(subdata.ddG, np.hstack([np.arange(-0.5, 5, 0.5), subdata.ddG.max()]),
                                                  precision=1, include_lowest=True)
            order = subdata.groupby('binned_ddG').first().index.tolist()
            num_annotations = pd.concat({name:group.annotation.value_counts() for name, group in subdata.groupby('binned_ddG')}).unstack().loc[order].fillna(0)
            fraction_annotation = (num_annotations.transpose()/num_annotations.sum(axis=1)).transpose()
            
            # add background expections
            num_background_annotations = data.loc[data.name.str.find('hPUM') == -1].annotation.value_counts().loc[annotation_order]
            expected_fractions = num_background_annotations/num_background_annotations.sum()
            fraction_annotation.loc['expected'] = expected_fractions

            # plot
            fraction_annotation.loc[:, annotation_order].plot(kind='bar', stacked=True, colors=annotation_colors, figsize=(3,3), width=0.6)
            plt.ylim(0, 1)
            plt.tight_layout()
            
            # plot enrichment relative to high ddG bin
            np.log2((fraction_annotation/fraction_annotation.loc[order[-1]]).loc[order, annotation_order]).plot(kind='bar', colors=annotation_colors, figsize=(3,3), width=0.8)
            
        
        elif args.mode == 'plot_versus_proximity_old':
            """Plot how the distribution of clip enrichment changes if there are nearby UGUA sites"""
            binedges = [0,  60,  150, 251]
            data.loc[:, 'binned_min_dist'] = pd.cut(data.min_bases_to_TGTA, binedges, include_lowest=True)
            
            subdata = data.loc[(data.tpm>0)&(data.annotation=="3' UTR")]
            subdata.loc[:, 'binned_ddG'] = pd.cut(subdata.ddG, np.hstack([np.arange(-0.5, 5, 0.5), subdata.ddG.max()]),
                                                  precision=1, include_lowest=True)
            order = subdata.groupby('binned_ddG').first().index.tolist()

            g = sns.factorplot(data=subdata, x='binned_ddG', hue='binned_min_dist', y='clip_signal_per_tpm', estimator=np.median,
                                     errwidth=1,capsize=1, linestyles='', marker='.', palette='viridis')
            plt.yscale('log'); plt.xticks(rotation=90); plt.subplots_adjust(bottom=0.35)
            ymax = subdata.groupby('binned_ddG').get_group(order[0]).clip_signal_per_tpm.median()            
            x = pd.Series({name:0.5*(group.ddG.max() + group.ddG.min()) for name, group in subdata.groupby('binned_ddG')})
            y = ymax*np.exp(x/seqmodel.get_ddG_conversion(temperature))
            
            plt.plot(np.arange(len(order)), y.loc[order], 'k--')
            ylim = [0.01, 20]
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
            
            
    elif args.mode == 'compare_to_pum12_kd' or args.mode == 'compare_to_pum2_oe':
        """Load supp data from NAR paper and compare sites."""
        
        # load the biomart ref
        biomart_data = pd.read_table('annotations/ensemble_gene_converter_biomart.txt', names=['gene_id', 'transcript_id', 'gene_name', 'refseq_id', 'refseq_nc'], header=0)
        biomart_data.loc[:, 'refseq_comb'] = [refseq_id if not str(refseq_id)=='nan' else refseq_nc for idx, refseq_id, refseq_nc in biomart_data.loc[:, ['refseq_id', 'refseq_nc']].itertuples()]
        if args.mode == 'compare_to_pum12_kd':
            col_name = 'gene_name'
        else:
            col_name = 'gene_id'
        data.loc[:, 'gene_name'] = pd.Series(biomart_data.groupby('refseq_comb').first().loc[data.refseq_id.dropna()][col_name].values, index=data.refseq_id.dropna().index)
                
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
            num_consensus_CDS = (group.loc[(group.annotation=="exon")].ddG < 0.5).sum()

            num_sites_2kc_3UTR = (group_3UTR.ddG < 2).sum()
            num_sites_between1and4kc_3UTR = ((group_3UTR.ddG < 4)&(group_3UTR.ddG >= 1)).sum()
            min_dG_3UTR = group_3UTR.ddG.min()
            occupancy_not3UTR = np.exp(-group.loc[group.annotation!="3' UTR"].ddG/RT).sum()
            occupancy_data[name] = pd.Series({'occupancy_3UTR':occupancy_3UTR,
                                              'occupancy_noflip_3UTR':occupancy_noflip_3UTR,
                                              'occupancy_not3UTR':occupancy_not3UTR,
                                              'min_ddG':min_dG_3UTR,
                                              'num_consensus_3UTR':num_consensus_3UTR,
                                              'num_consensus_CDS':num_consensus_CDS,
                                              'num_sites_2kc_3UTR':num_sites_2kc_3UTR,
                                              'num_sites_between1and4kc_3UTR':num_sites_between1and4kc_3UTR})
        occupancy_data = pd.concat(occupancy_data).unstack()
        occupancy_data.loc[occupancy_data.min_ddG.isnull(), 'min_ddG'] = 10 #kcal.mol
        
        # load expression data
        if args.mode == 'compare_to_pum12_kd':
            expression_data = pd.read_csv('annotations/nar-01280/supp_table4.csv')
            expression_data.loc[:, 'lfc'] = expression_data.lfc.replace('#NAME?', np.nan).astype(float).replace(np.inf, np.nan)        
            expression_data.loc[:, 'log_fpkm_cntrl'] = np.log10(expression_data.FPKM_NTC)
            expression_fpkm_bins = [-np.inf, 1, 1.5]
        elif args.mode == 'compare_to_pum2_oe':
            expression_data = pd.read_table('annotations/GSE75440/GSe75440_PUM2edgeR.txt.gz')
            expression_data.rename(columns={col:col.lower().replace(' ', '_').replace('.', '_')  for col in expression_data}, inplace=True)  
            expression_data.loc[:, 'lfc'] = expression_data['logfc_pum2/gfp']
            expression_data.loc[:, 'log_fpkm_cntrl'] = np.log10(expression_data.loc[:, 'gfp_#1_fpkm':'gfp_#3_fpkm']+0.01).mean(axis=1)
            #expression_data.loc[:, 'gene_id'] = expression_data.gene
            #expression_data.loc[:, 'gene'] = expression_data.genename
            
            expression_data.loc[:, 'sig_up'] = (expression_data.adj_pval_tgw < 1E-2)&(expression_data.lfc > 0)
            expression_data.loc[:, 'sig_down'] = (expression_data.adj_pval_tgw < 1E-2)&(expression_data.lfc < 0)
            #expression_data.set_index('gene', inplace=True)
            expression_fpkm_bins = [-np.inf, 1, 3]
   
        # find roc
        roc_curves = {}
        for name in occupancy_data:
            
            expression_data.loc[:, 'occupancy'] = occupancy_data.loc[expression_data.gene, name].values
            expression_data_sub = expression_data.dropna(subset=['lfc', 'occupancy']).copy()
            if name=='min_ddG':
                predicted_up=False
            else:
                predicted_up = True
            roc_curves[('any', name)] =  processing.find_roc_data(expression_data_sub, 'occupancy', 'sig_down')
            #for expression_threshold in expression_fpkm_bins:
            #    roc_curves[(expression_threshold, name)] =  processing.find_roc_data(expression_data_sub.loc[expression_data_sub.log_fpkm_cntrl >= expression_threshold], 'occupancy', 'sig_up', predicted_up=predicted_up, )
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
        
 