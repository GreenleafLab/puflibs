from fittinglibs import processresults, fileio


annotatedClusters = fileio.loadFile('subset_variants.CPannot.pkl')
bindingSeries = fileio.loadFile('sarah_try0/AM4FH_ALL_Bottom_filtered_reduced_normalized.CPseries.pkl')
concentrations = fileio.loadFile('concentrations.txt')
fmaxDist = fileio.loadFile('sarah_try0/AM4FH_ALL_Bottom_filtered_reduced_normalized.fmaxdist.p')
fitResults = fileio.loadFile('sarah_try0/AM4FH_ALL_Bottom_filtered_reduced_normalized.CPfitted.pkl')

filenameDict = {'original': 'AM4FH_ALL_Bottom_filtered_reduced_normalized.CPvariant',
'no weights': 'AM4FH_ALL_Bottom_filtered_reduced_normalized_noweights.CPvariant',
'fmin float': 'AM4FH_ALL_Bottom_filtered_reduced_normalized_fminfloat.CPvariant',
'no weights, slope float': 'AM4FH_ALL_Bottom_filtered_reduced_normalized_noweights_slopeflaot.CPvariant',
'fmin float, slope float': 'AM4FH_ALL_Bottom_filtered_reduced_normalized_fminfloat_slopeflaot.CPvariant',
'fmax float, fmin float, no weights': 'AM4FH_ALL_Bottom_filtered_reduced_normalized_fmaxfloat_fminfloat_noweights.CPvariant',
'fmax float, fmin float, no weights, slope float': 'AM4FH_ALL_Bottom_filtered_reduced_normalized_fmaxfloat_fminfloat_noweights_slopefloat.CPvariant'}

accepted_variants = np.arange(20138.)
order = ['original', 'no weights', 'fmin float', 'no weights, slope float', 'fmin float, slope float', 'fmax float, fmin float, no weights', 'fmax float, fmin float, no weights, slope float']

affinityData = {}
names = []
for name, filename in filenameDict.items():
    names.append(name)
    variant_table = fileio.loadFile(os.path.join('sarah_try0/', filename))
    affinityData[name] = (processresults.perVariant(variant_table.loc[accepted_variants].dropna(how='all'), annotatedClusters, bindingSeries, cluster_table=fitResults, x=concentrations))

variant_table = fileio.loadFile('sarah_try1/AM4FH_ALL_Bottom_filtered_reduced_normalized.CPvariant')
affinityData['no weights, nonlinear'] = (processresults.perVariant(variant_table.loc[accepted_variants].dropna(how='all'), annotatedClusters, bindingSeries, cluster_table=fitResults, x=concentrations))

variant_table = fileio.loadFile('sarah_try1/AM4FH_ALL_Bottom_filtered_reduced_normalized_dGns_fixed.CPvariant')
affinityData['no weights, nonlinear fixed'] = (processresults.perVariant(variant_table.loc[accepted_variants].dropna(how='all'), annotatedClusters, bindingSeries, cluster_table=fitResults, x=concentrations))
  
expt = processresults.manyFlows(affinityData, names)
expt.plotAllInitVsFinal(order=order, vmin=0.29, vmax=1.29)

### compare two sets
dGmat = pd.concat([affinityData['no weights'].variant_table.dG.rename('slope_fixed'), affinityData['no weights, slope float'].variant_table.dG.rename('slope_float')], axis=1).dropna(how='all')
fmaxfixedMat = pd.concat([affinityData['no weights'].variant_table.flag.rename('slope_fixed'), affinityData['no weights, slope float'].variant_table.flag.rename('slope_float')], axis=1).loc[dGmat.index]
slopeMat = pd.concat([pd.Series(6E-3, index=dGmat.index).rename('slope_fixed'), affinityData['no weights, slope float'].variant_table.slope.rename('slope_float')], axis=1).loc[dGmat.index]

plt.figure(); plt.hexbin(dGmat.slope_fixed, dGmat.slope_float, cmap='copper', extent=[-12.5, -7, -12.5, -7], gridsize=200, mincnt=1)
plt.figure(); plt.scatter(dGmat.slope_fixed, dGmat.slope_float, c=fmaxfixedMat.slope_fixed, cmap='bwr')

# plot scatterplot, colored by fit slope
plt.figure(figsize=(4,3)); im=plt.scatter(mat.slope_fixed, mat.slope_float, c=slopeMat.slope_float, cmap='viridis', vmin=0, vmax=0.0012, s=10); sns.despine(); fix_axes(plt.gca()); plt.colorbar(im); plt.xlim(-12.5, -8); plt.ylim(-12.5, -8); plt.xticks(np.arange(-12, -7)); plt.yticks(np.arange(-12, -7)); plt.ylabel('slope float'); plt.xlabel('slope fixed'); plt.tight_layout()


### compare two sets
dGmat = pd.concat([affinityData['no weights'].variant_table.dG.rename('slope_fixed'), affinityData['no weights, nonlinear'].variant_table.dG.rename('nonlinear')], axis=1).dropna(how='all')
fmaxfixedMat = pd.concat([affinityData['no weights'].variant_table.flag.rename('slope_fixed'), affinityData['no weights, nonlinear'].variant_table.flag.rename('nonlinear')], axis=1).loc[dGmat.index]
dGns = affinityData['no weights, nonlinear'].variant_table.dG_ns.rename('nonlinear').loc[dGmat.index]

# plot scatterplot, colored by fit slope
plt.figure(figsize=(4,3)); im=plt.scatter(dGmat.slope_fixed, dGmat.nonlinear, c=dGns, cmap='viridis', s=10, vmin=-8, vmax=60); sns.despine(); fix_axes(plt.gca()); plt.colorbar(im); plt.xlim(-12.5, -8); plt.ylim(-12.5, -8); plt.xticks(np.arange(-12, -7)); plt.yticks(np.arange(-12, -7)); plt.ylabel('slope float'); plt.xlabel('slope fixed'); plt.tight_layout()

### compare two sets
dGmat = pd.concat([affinityData['no weights'].variant_table.dG.rename('slope_fixed'), affinityData['no weights, nonlinear fixed'].variant_table.dG.rename('nonlinear')], axis=1).dropna(how='all')

# plot scatterplot, colored by fit slope
plt.figure(figsize=(4,3)); im=plt.scatter(dGmat.slope_fixed, dGmat.nonlinear, c=dGns, cmap='viridis', s=10, vmin=-8, vmax=60); sns.despine(); fix_axes(plt.gca()); plt.colorbar(im); plt.xlim(-12.5, -8); plt.ylim(-12.5, -8); plt.xticks(np.arange(-12, -7)); plt.yticks(np.arange(-12, -7)); plt.ylabel('slope float'); plt.xlabel('slope fixed'); plt.tight_layout()

### compare two sets
dGmat = pd.concat([affinityData['no weights, slope float'].variant_table.dG.rename('slope_float'), affinityData['no weights, nonlinear fixed'].variant_table.dG.rename('nonlinear_fixed')], axis=1).dropna(how='all')
slopeMat = affinityData['no weights, slope float'].variant_table.slope.rename('slope_float').loc[dGmat.index]

dGns = affinityData['no weights, nonlinear'].variant_table.dG_ns.rename('nonlinear').loc[dGmat.index]

# plot scatterplot, colored by fit slope
plt.figure(figsize=(4,3)); im=plt.scatter(dGmat.slope_float, dGmat.nonlinear_fixed, c=slopeMat, cmap='viridis',  vmin=0, vmax=0.0012, s=10); sns.despine(); fix_axes(plt.gca()); plt.colorbar(im); plt.xlim(-12.5, -8); plt.ylim(-12.5, -8); plt.xticks(np.arange(-12, -7)); plt.yticks(np.arange(-12, -7)); plt.ylabel('slope float'); plt.xlabel('slope fixed'); plt.tight_layout()
