
for out_dir in "${out_dirs[@]}"
	do
	  	a=$out_dir; prefix="$a"

		let i=0		
		for f in ` ls $out_dir/*.dat ` 
  			do
	  			data_files[$i]=` basename $f `
				i=$(($i+1))
		done

		cd $out_dir
			source $this_script_dir/gp_loop_over_dat.sh
		cd $this_script_dir
done

