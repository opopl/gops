
		#print header  {{{
		lst="------------------------------------------------------------"

		l[0]="XYZ file: $file0; Frame: $i_frame; Output image file: $pof" 
		l[1]="Parameters: E= $energy; F= $force; N=$nst; T=$temp; R=$radius" 
		l[2]="Rotation angles: X= $x_angle; Y= $y_angle; Z= $z_angle" 
		l[3]="Zoom: $zoom" 

		for line in "${l[@]}"
			do
				printf "%s\n" "$line"
	      	done	
		
		case "$print_text" in
		  	"all") 	
				jtext="" 
				for line in "${l[@]}"
					do
						jtext=` printf "%s  " "$jtext" "$line" "\n"  `
				done
				;;
			"no") jtext="" ;;
			"ra") jtext="$x_angle $y_angle $z_angle" ;;
		esac
	  #}}}
		
		pic="$pics_dir/pic.$pic_ext"
	   	source "$this_script"_j.sh 		

		jmol -ions $r/t.jmol >& jmol.log

		if [ $angle_count -eq 0 ]; then
		  	cp $pic $output_pic_file 
		else
			gifsicle $pic $output_pic_file > n 
			mv n $output_pic_file 
			rm -f $pic
		fi
		rm -f $pic

		pic_file_size=$(stat -c%s "$output_pic_file")
		
		echo "Total image file size: $pic_file_size"

