cd /lab/sarah/pufProject/AMEMH/analysis/binding_curves
annot="../amemh_all_tiles_reduced.CPannot.pkl"
variant_table="25C_refit/AMEMH_ALL_Bottom_filtered_reduced_normalized.init.CPvariant.pkl"
binding_series="25C_refit/AMEMH_ALL_Bottom_filtered_reduced_normalized.CPseries.pkl"
fmaxdist="25C_refit/AMEMH_ALL_Bottom_filtered_reduced_normalized.fmaxdist.p"

out=$(basename $variant_table .init.CPvariant.pkl)""
python ~/array_fitting_tools/bin/bootStrapFits.py --no_weights -v $variant_table -a $annot -b $binding_series -c concentrations.txt -f $fmaxdist --variants variants.txt -out sarah_try1/output.no_weights