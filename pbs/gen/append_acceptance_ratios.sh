
#grep $dn/GMIN_out 'Acceptance ratio' | awk '{ print $7, print $9}' >> $output_dir/ac_ratios 

# how to include info about STEP as well ???
