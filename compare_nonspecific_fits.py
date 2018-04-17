
# load the different fits
fitParams = fileio.loadFile('25C_refit2_subset/AMEMH_ALL_Bottom_filtered_reduced_normalized_subset_median.fitParameters.p') # binding curve nonlinear
fitParamsOG = fileio.loadFile('25C_refit2_subset/AMEMH_ALL_Bottom_filtered_reduced_normalized_subset_median_no_dGns.fitParameters.p') # binding curve

# and results
binding_series = fileio.loadFile('25C_refit2_subset/AMEMH_ALL_Bottom_filtered_reduced_normalized_subset_median.CPseries.pkl')
fitted_vals = fileio.loadFile('25C_refit2_subset/AMEMH_ALL_Bottom_filtered_reduced_normalized_subset_median.CPfitted.pkl')
fitted_vals_fmax = fileio.loadFile('25C_refit2_subset/AMEMH_ALL_Bottom_filtered_reduced_normalized_subset_median_fmax_fixed_1.18.CPfitted.pkl')
fitted_vals_fixed = fileio.loadFile('25C_refit2_subset/AMEMH_ALL_Bottom_filtered_reduced_normalized_subset_median_dGns_fixed_-8.69.CPfitted.pkl')
fitted_vals_og = fileio.loadFile('25C_refit2_subset/AMEMH_ALL_Bottom_filtered_reduced_normalized_subset_median_no_dGns.CPfitted.pkl')
variant_table = fileio.loadFile('25C_refit2/AMEMH_ALL_Bottom_filtered_reduced_normalized_dGns_fixed_-8.69_new.CPvariant')
# compare
compare_results = {}
for param_name in ['dG', 'fmax', 'fmin', 'rmse']:
    compare_results[param_name] = pd.concat([data.loc[:, param_name].rename(name) for name, data in zip(['float', 'fixed', 'none', 'fmax_fixed'], [fitted_vals, fitted_vals_fixed, fitted_vals_og, fitted_vals_fmax])], axis=1)
compare_results = pd.concat(compare_results, axis=1)

# plot
xlim = [-14, -6]
plt.figure(figsize=(3,3))
plt.scatter(compare_results.loc[:, ('dG', 'none')], compare_results.loc[:, ('dG', 'float')], c=compare_results.loc[:, ('rmse', 'float')], cmap='coolwarm', vmin=0, vmax=0.7, marker='.')
plt.xlim(xlim)
plt.ylim(xlim)
plt.plot(xlim, xlim, 'k-', linewidth=1)
plt.xlabel('no nonspecific term (kcal/mol)')
plt.ylabel('nonspecific term float (kcal/mol)')
plt.tight_layout()

plt.figure(figsize=(3,3))
plt.scatter(compare_results.loc[:, ('dG', 'none')], compare_results.loc[:, ('dG', 'fixed')], c=compare_results.loc[:, ('rmse', 'fixed')], cmap='coolwarm', vmin=0, vmax=0.7, marker='.')
plt.xlim(xlim)
plt.ylim(xlim)
plt.plot(xlim, xlim, 'k-', linewidth=1)
plt.xlabel('no nonspecific term (kcal/mol)')
plt.ylabel('nonspecific term fixed (kcal/mol)')
plt.tight_layout()

plt.figure(figsize=(3,3))
plt.scatter(compare_results.loc[:, ('dG', 'fixed')], compare_results.loc[:, ('dG', 'float')], c=compare_results.loc[:, ('rmse', 'float')], cmap='coolwarm', vmin=0, vmax=0.7, marker='.')
plt.xlim(xlim)
plt.ylim(xlim)
plt.plot(xlim, xlim, 'k-', linewidth=1)
plt.xlabel('nonspecific term fixed (kcal/mol)')
plt.ylabel('nonspecific term float (kcal/mol)')
plt.tight_layout()


plt.figure(figsize=(3,3))
plt.scatter(compare_results.loc[:, ('dG', 'fixed')], compare_results.loc[:, ('dG', 'fmax_fixed')], c=compare_results.loc[:, ('rmse', 'fmax_fixed')], cmap='coolwarm', vmin=0, vmax=0.7, marker='.')
plt.xlim(xlim)
plt.ylim(xlim)
plt.plot(xlim, xlim, 'k-', linewidth=1)
plt.xlabel('nonspecific term fixed (kcal/mol)')
plt.ylabel('nonspecific term float (kcal/mol)')
plt.tight_layout()

plt.figure(figsize=(3,3))
plt.scatter(compare_results.loc[:, ('dG', 'float')], compare_results.loc[:, ('dG', 'fmax_fixed')], c=compare_results.loc[:, ('rmse', 'fmax_fixed')], cmap='coolwarm', vmin=0, vmax=0.7, marker='.')
plt.xlim(xlim)
plt.ylim(xlim)
plt.plot(xlim, xlim, 'k-', linewidth=1)
plt.xlabel('nonspecific term fixed (kcal/mol)')
plt.ylabel('nonspecific term float (kcal/mol)')
plt.tight_layout()

# plot the nonlinear term
plt.figure(figsize=(3,3))
plt.hexbin(fitted_vals.loc[:, 'dG'], fitted_vals.loc[:, 'dGns'], cmap='copper', extent=[-14, -6, -12, -4], mincnt=1, bins='log') 


# plot the nonlinear term
plt.figure(figsize=(3,3))
plt.hexbin(fitted_vals_fmax.loc[:, 'dG'], fitted_vals_fmax.loc[:, 'dGns'], cmap='copper', extent=[-14, -6, -12, -4], mincnt=1, bins='log') 


# plot the nonlinear term
plt.figure(figsize=(3,3))
plt.hexbin(variant_table.loc[:, 'dG_init'], variant_table.loc[:, 'dGns_init'],  cmap='copper', extent=[-14, -6, -12, -4], mincnt=1, bins='log') 
plt.figure(figsize=(3,3))
plt.hexbin(variant_table.loc[fitted_vals.index, 'dG_init'], variant_table.loc[fitted_vals.index, 'dGns_init'],  cmap='copper', extent=[-14, -6, -12, -4], mincnt=1, bins='log')

plt.figure(figsize=(3,3))
plt.scatter(variant_table.loc[fitted_vals.index, 'dG_init'], variant_table.loc[fitted_vals.index, 'dGns_init'],marker='.')
plt.figure(figsize=(3,3))
plt.scatter(fitted_vals.loc[:, 'dG'], fitted_vals.loc[:, 'dGns'], marker='.')
plt.figure(figsize=(3,3)); tectplots.distplot(variant_table.loc[fitted_vals.index, 'dGns_init']); tectplots.distplot(variant_table.dGns_init); plt.axvline(-8.69)

# plot the difference between the two fits

diff_vec = (compare_results.loc[:, ('dG', 'float')] - compare_results.loc[:, ('dG', 'fixed')])
index = (compare_results.loc[:, ('rmse', 'float')]<0.7)&(compare_results.loc[:, ('dG', 'float')]<-11)
index2 = (compare_results.loc[:, ('dG', 'fixed')]<-12)&(compare_results.loc[:, ('dG', 'fixed')]>-12.1)
plt.figure(figsize=(3,3))
tectplots.distplot(diff_vec.loc[index]); tectplots.distplot(diff_vec.loc[index&index2])

binned_vals = pd.Series(np.digitize(diff_vec.loc[index&index2], [-7, -0.35, -0.1, 0.1, 0.35, 7] ), index=diff_vec.loc[index&index2].index)
for i in binned_vals.unique():
    idxs = binned_vals.loc[binned_vals==i].index.tolist()
    for idx in np.random.choice(idxs, size=min(5, len(idxs)), replace=False):
        plt.figure(figsize=(3,3))
        fitParamsOG.plot_fit(results=fitted_vals_og.loc[idx], y=binding_series.loc[idx], label='no nonspecific term', marker=',')
        fitParams.plot_fit(results=fitted_vals_fixed.loc[idx], y=binding_series.loc[idx], label='nonspecific term fixed', marker=',')
        fitParams.plot_fit(results=fitted_vals.loc[idx], y=binding_series.loc[idx], label='nonspecific term float', marker=',')
        fitParams.plot_fit( y=binding_series.loc[idx], color='k',marker='.')
        plt.xscale('log')
        plt.xlabel('concentration (nM)')
        plt.ylabel('fluorescence')
        #plt.legend()
        plt.tight_layout()
        plt.title('%d; diff=%4.2f'%(idx, diff_vec.loc[idx]))
        plt.savefig('25C_refit2_subset/figures/binding_curve.%d.binned_diff_%d.dG_-12.pdf'%(idx, i))
