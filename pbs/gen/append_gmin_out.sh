
export outf=$output_dir/gmin_out
eo_s "$xyz_comment_string" 
cat $dn/GMIN_out >> $outf
#cat $dn/GMIN_out | grep 'WARNING'
