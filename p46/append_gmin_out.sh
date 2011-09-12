
export outf=$output_dir/gmin_out
echo "==========================" >> $outf
eo_s "$var_full_string" 
echo "==========================" >> $outf
cat $dn/GMIN_out >> $outf


