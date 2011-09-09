
fouts=( energy.dat markov.dat ) 

for fout in "${fouts[@]}"
	do
		echo "# @@ $xyz_comment_string" >> $output_dir/$fout
		cat $dn/$fout >> $output_dir/$fout
done
