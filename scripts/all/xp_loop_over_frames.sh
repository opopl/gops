#!/bin/bash - 

# loop over frames {{{

for ((i_frame=$frame_min; i_frame<=$frame_max; i_frame++ )); do
	  	# j.sh produces t.jmol
		echo "File: $file Frame index: $i_frame"
		echo "====================================="
		xyz_header=${xyz_headers[$i_frame]}
		
		mod="$file0"_"$i_frame"
		xyzfile="$tmp_dir/$mod.xyz"

		gawk 'NR==start, NR==end' start=$((48*$i_frame-47)) end=$((48*$i_frame)) $full_file >& $xyzfile 

		echo "XYZ file $file0> XYZ frame header: $xyz_header"
		echo "====================================="

		energy=` echo $xyz_header | awk '{ print $2 }' | sed 's/;$//' ` 
		force=` echo $xyz_header | awk '{ print $4 }' | sed 's/;$// ' ` 
		nsteps=` echo $xyz_header | awk '{ print $6 }' | sed 's/;$//'  ` 
	
		if [  $reset_pof -ne 1 ]; then	
			pof="f"$force"_1"."$pic_ext"
	        fi

		pofo=$pof

		let i=1

		if [ $append_force -eq 1 ]; then 
			while [ -f $pof ]; do 
	  			pof=${pofo%.1.$pic_ext}_$i.$pic_ext
				i=$(($i+1))
			done
	        fi
		
		output_pic_file="$pics_dir/$pof"
		#touch output_pic_file
		
		case "$nsteps" in
		  	1000000) nst="1E6" ;;
		  	100000) nst="1E5" ;;
		  	10000) nst="1E4" ;;
		  	1000) nst="1E3" ;;
		esac

		radius=` echo $xyz_header | awk '{ print $12 }' | sed 's/;$//' | sed 's/^R=//'  ` 
		temp=` echo $xyz_header | awk '{ print $11 }' | sed 's/;$//' | sed 's/^T=//'  ` 
		
		output_pic_file_old=$output_pic_file

		if [ $use_opt -eq 1 ]; then  
			use use_opt
	      	fi
		
		use rotate
	  done
	#}}}

