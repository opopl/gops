
export outf=$output_dir/pdf

eo_s "$xyz_comment_string" 

cat $dn/pdf >> $outf

phrase="mu/sig of t, Q, V  are"

# recall, here f is the force 
#cat $dn/pdf  | sed "s/^[a-z/A-Z, ]*/$f  /" >> $output_dir/pdff.dat
#cat $dn/pdf  | awk '{ print $7}' >> $output_dir/pdff.dat


cat $dn/pdf | awk ' {
for (i=7; i<=NF; i++)
  printf("%s ", $i)
  printf("\n") # CR at end of line
} ' | sed "s/^/$var_string /">> $output_dir/pdff.dat
