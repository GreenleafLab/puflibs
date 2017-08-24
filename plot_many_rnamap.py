#!/usr/bin/env python
""" Make figures for paper.

Sarah Denny """

##### IMPORT #####
import numpy as np
import pandas as pd
import os
import argparse
import pickle
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import scipy.stats as st
import scipy.cluster.hierarchy as sch
from sklearn import metrics
from scikits.bootstrap import bootstrap
from scipy.cluster.vq import kmeans2, vq
from sklearn.manifold import TSNE
import scipy.spatial.distance as ssd
import lmfit
import itertools
import functools
from statsmodels.sandbox.stats.multicomp import multipletests
import hjh.mutations
from fittinglibs import processresults, fileio, seqfun, fitting, objfunctions
from fittinglibs.plotting import fix_axes
from tectolibs import (resultscompare, exptplots, tectplots, clustering)

RT = 0.593 # at 25 deg C

### MAIN ###
#set up command line argument parser
parser = argparse.ArgumentParser(description='bootstrap on off rate fits')
parser.add_argument('--data', default="data_refit/analysis/filename_table.dat",
                   help='file giving filenames')
parser.add_argument('--key', nargs='*',
                   help='key of experiment to look at')
parser.add_argument('--mode', 
                   help='what to do with it')
parser.add_argument('--out', 
                   help='indication of output file')
parser.add_argument('-f', '--force_overwrite', action="store_true", 
                   help='whether to overwrite output')
args = parser.parse_args()
filename_table = pd.read_table(args.data, index_col=0)

dGmax = -9.913941 # 

if args.mode == 'find_results_table':
    """Using the variant table, find the results table."""
    if args.out is None:
        print 'Error: please supply --out parameter. Something that looks like: '
        print 'paper/01_expt/results_tables/flow_3455.151204.error_scaled.results.pkl'
        sys.exit()
        
    if os.path.exists(args.out) and not args.force_overwrite:
        print 'Error: filename %s already exists. Use -f to overwrite.'%args.out
        sys.exit()
    
    key = args.key[0]
    if not os.path.exists(filename_table.loc[key, 'variant_table']):
        print 'Error: need to supply variant table filename for key %s in filename table %s.'%(key, args.data)
        sys.exit()

    variant_table = fileio.loadFile(filename_table.loc[key, 'variant_table'])
    if 'numTests' not in variant_table.columns.tolist():
        variant_table.rename(columns={'num_tests':'numTests'}, inplace=True)
    affinity_data = exptplots.PerVariant(variant_table=variant_table)
    result_table = affinity_data.getResultsFromVariantTable()
    fileio.saveFile(args.out, result_table)

if args.mode == 'combine_results_table':
    """Combine two variant_tables to form a results table."""
    if args.out is None:
        print 'Error: please supply --out parameter. Something that looks like: '
        print 'paper/01_expt/results_tables/flow_3455.151204.error_scaled.results.pkl'
        sys.exit()
        
    if os.path.exists(args.out) and not args.force_overwrite:
        print 'Error: filename %s already exists. Use -f to overwrite.'%args.out
        sys.exit()
    
    keys = args.key
    for key in keys:
        if not os.path.exists(filename_table.loc[key, 'variant_table']):
            print 'Error: need to supply variant table filename for key %s in filename table %s.'%(key, args.data)
            sys.exit()

    variant_tables = [fileio.loadFile(filename_table.loc[key, 'variant_table']) for key in keys]
    for variant_table in variant_tables:
        if 'numTests' not in variant_table.columns.tolist():
            variant_table.rename(columns={'num_tests':'numTests'}, inplace=True)

    result_table = processresults.getResultsFromVariantTables(variant_tables, offset=0)
    fileio.saveFile(args.out, result_table)
    
    
if args.mode == 'scale_variant_table_error':
    """Scale the error bars of variant table by known funciton, and save as variant table."""
    if args.out is None:
        print 'Error: please supply --out parameter. Something that looks like: '
        print 'paper/01_expt/variant_tables/flow_WC0.141111.error_scaled.CPvariant'
        sys.exit()
        
    if os.path.exists(args.out) and not args.force_overwrite:
        print 'Error: filename %s already exists. Use -f to overwrite.'%args.out
        sys.exit()
    
    key = args.key[0]
    if not os.path.exists(filename_table.loc[key, 'variant_table']):
        print 'Error: need to supply variant table filename for key %s in filename table %s.'%(key, args.data)
        sys.exit()
        
    variant_table = fileio.loadFile(filename_table.loc[key, 'variant_table'])
    affinity_data = exptplots.PerVariant(variant_table=variant_table)
    new_variant_table = affinity_data.correct_variant_table(amplitude=0.009860, exponent=-1.748300, c=0.929585)
    fileio.saveFile(args.out, new_variant_table)
    """
    results_table = affinityData.getResultsFromVariantTable()
    results = resultscompare.SingleResult(results_table)
    results_table_corr = results.correct_results()
    """

elif args.mode == 'plot_error_figs':
    """Compare the two replicate datasets."""
    lib_char = fileio.loadFile(filename_table.loc['libchar1'].libchar)

    names = ['hPUM2_1_uncorr', 'hPUM2_2_uncorr', 'hPUM2_1', 'hPUM2_2']
    results_dict = {key:fileio.loadFile(filename_table.loc[key].results_table) for key in names}
    compare = resultscompare.AllResults(results_dict, lib_char, names)
    # divide into intervals
    name_dict = {'uncorr':['hPUM2_1_uncorr', 'hPUM2_2_uncorr'], 'corr':['hPUM2_1', 'hPUM2_2']}
    indices = compare.find_indices(mode='err_bins', names=name_dict['uncorr'], cutoff_mean=True, dGmax=dGmax)
    index = list(itertools.chain(*indices.values()))
    # find fdr in different intervals
    #fdr_dict = compare.find_frac_passing_fdr(mode='dGmax', names=names, return_fdr=True)
    for key, names in name_dict.items():
        # calculate zscores of different between groups
        offset = compare.find_offsets(mode='dGmax', names=names, dGmax=dGmax)[0]['dGmax']
        #offset, index = offsets['dGmax'], indices['dGmax']
        sigma_both = compare.find_expt_sigma(names=names).rename('sigma')
        zscores = ((compare.find_expt_diff(names=names) - offset)/sigma_both).rename('zscores')
        
        # plot zscores
        x = np.linspace(-10, 10, 100)
        tectplots.figure(); sns.distplot(zscores.loc[index].dropna()); plt.plot(x, st.norm.pdf(x), 'k--')
        
        # find ratios
        mat = tectplots.reshape_data_by_indices(indices, zscores).set_index(['index_key', 'index'])
        fold_sigmas = {}
        for name, group in mat.groupby(level=0):
            params, stderr = tectplots.fitGaussian(group.zscores)
            lb, x, ub = sigma_both.loc[indices[name]].dropna().quantile([0, 0.5, 1])
            fold_sigmas[name] = pd.Series([lb, ub, x, params.sigma, stderr.sigma, params.center],
                index=['lb', 'ub', 'x', 'y', 'y_stderr', 'offset'])
        fold_sigmas = pd.concat(fold_sigmas).unstack()
        
        # fit to exponential
        fitParameters = pd.concat({'amplitude':pd.Series({'initial':0.7}), 'exponent':pd.Series({'initial':-0.25}), 'c':pd.Series({'initial':0})}).unstack(level=0).fillna(True)
        params = fitting.convertFitParametersToParams(fitParameters)
        results = fitting.minimize(objfunctions.powerlaw, params, args=(fold_sigmas.dropna().x,), kws={'y':fold_sigmas.dropna().y})
        final_params = fitting.returnResultsFromParams(params, results, fold_sigmas.y)
        
        # plot
        plt.figure(figsize=(3,3)); plt.scatter(fold_sigmas.x, fold_sigmas.y, c='0.5', edgecolor='w');
        more_x = np.linspace(0, 3)
        plt.plot(more_x, objfunctions.powerlaw(params, more_x), 'r', linewidth=1)
        plt.xlim(0, 2.5)
        plt.ylim(0, 2)
        plt.axhline(1, linewidth=1, linestyle=':', color='0.5')
        sns.despine()
    
    

if args.mode == 'compare_replicates':
    """copmare the two replicat datasets and find FDR etc."""
    lib_char = fileio.loadFile(filename_table.loc['libchar1'].libchar)

    names = ['hPUM2_1', 'hPUM2_2']
    results_dict = {key:fileio.loadFile(filename_table.loc[key].results_table) for key in names}
    compare = resultscompare.AllResults(results_dict, lib_char, names)
    
    offset, index = [val['dGmax'] for val in compare.find_offsets(mode='dGmax', names=names, dGmax=dGmax, error_cutoff=True)]
    tectplots.figure(); plt.hexbin(compare.param_values.loc[index, names[0]], compare.param_values.loc[index, names[1]], cmap='copper', mincnt=1)
    xlim = np.array([-14, -9.9])
    plt.plot(xlim, xlim+offset)

if args.mode == 'plot_single_muts':
    """compare the single mutants"""
    lib_char = fileio.loadFile(filename_table.loc['libchar1'].libchar)
    single_muts = pd.read_csv(filename_table.loc['singlemuts', 'data'], index_col=0)
    single_muts.index = single_muts.unique_index
    single_muts.dropna(subset=['unique_index'], inplace=True)
    single_muts.sort_index(inplace=True)
    # process
    single_muts.loc[:, 'ss_25'] = single_muts.loc[:, '25_intr_dG'] - single_muts.loc[:, '25_dG']
    single_muts.loc[:, 'ss_37'] = single_muts.loc[:, '37_intr_dG'] - single_muts.loc[:, '37_dG']
    single_muts.drop(single_muts.loc[:, '25_dG':'37_num_tests'].columns.tolist(), axis=1, inplace=True)

    # laod data    
    names = ['hPUM2', 'hPUM2_37']
    results_dict = {key:fileio.loadFile(filename_table.loc[key].results_table) for key in names}
    compare = resultscompare.AllResults(results_dict, lib_char, names)
    
    # add data
    for col in names:
        single_muts.loc[:,col] = compare.param_values.loc[single_muts.index, col]
    
    # find ddG from consensus
    consensus = single_muts.set_index(['VarPosition', 'VarSeq']).loc[(1,'U'), names].mean()
    for col in names:
        single_muts.loc[:,col + '_rel'] = compare.param_values.loc[single_muts.unique_index, col].values - consensus.loc[col]
    
    # sort by mut/pos 
    indices_series = single_muts.set_index(['VarPosition', 'VarSeq', 'PUF_SCAFFOLD']).unique_index.dropna()
    indices = {name:list(group.values) for name, group in indices_series.groupby(level=[0,1])}
    
    # find sig diff from scaffold mean for each mut
    param_errors = compare.param_errors.loc[single_muts.index]
    param_values = compare.param_values.loc[single_muts.index]
    for col, temp in zip(names, ['25', '37']):
        param_errors.loc[:, col+'_ss'] = compare.param_errors.loc[single_muts.index, col]
        param_values.loc[:, col+'_ss'] = compare.param_values.loc[single_muts.index, col] + single_muts.loc[:, 'ss_'+temp]

    zscores_all, diff_all, sigma_all = tectplots.get_zscores_from_mean(param_values, param_errors, indices)

    # plot those called as sig siff
    fdr = pd.concat({col:tectplots.get_fdr(zscores_all.loc[:, col]) for col in zscores_all}, axis=1)
    sig_diff = fdr < 0.1
    for col in names:
        single_muts.loc[:, col + '_sig'] = sig_diff.reset_index(level=0).loc[single_muts.unique_index, col].values
    for col in names:
        sns.factorplot(data=single_muts, col='VarPosition', x='VarSeq', hue='PUF_SCAFFOLD', y=col + '_rel',
                       kind='bar', order=['A', 'C', 'G', 'U'], size=2.5, aspect=0.8,
                       hue_order=['S1a', 'S1b', 'S2a', 'S2b']); plt.ylim(-0.5, 4.5)
        sns.factorplot(data=single_muts.loc[single_muts.loc[:, col+'_sig'].fillna(False)], col='VarPosition', x='VarSeq', hue='PUF_SCAFFOLD', y=col+'_rel',
                       kind='bar', order=['A', 'C', 'G', 'U'], size=2.5, aspect=0.8,
                       hue_order=['S1a', 'S1b', 'S2a', 'S2b']); plt.ylim(-0.5, 4.5)
        
    # plot difference from expected
    for col in [x+'_ss' for x in names]:
        sim_values = [st.norm.rvs(scale=i) for  i in np.random.choice(sigma_all.loc[:, col].dropna(), size=1000)]
        tectplots.figure()
        tectplots.plot_two_dists(diff_all.loc[:, col].dropna(), sim_values)
        plt.xlim(-2, 2)
        plt.ylim(0, 3.6)
        plt.tight_layout()
    
    # correct for structure
    param_errors = compare.param_errors
    param_values = compare.param_values.loc[single_muts.unique_index] + pd.concat({'hPUM2':single_muts.loc[:, '25_intr_dG'] - single_muts.loc[:, '25_dG'] ,
                                                                                   'hPUM2_37':single_muts.loc[:, '37_intr_dG'] - single_muts.loc[:, '37_dG']}, axis=1)
        