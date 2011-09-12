
f=$xyz_dir/$1.xyz
	
rgyr_file=$rgyr_dir/rgyr_$1.dat
echo "#" >& $rgyr_file

eh "Input XYZ file: $1.xyz"

fo=$f.tmp
cp $f $fo

sz_fo=` wc -c $fo | awk '{ print $1 }' ` 

xyz=tmp342
this_frame=tmp3424

let i_frame=1

# loop over all the frames in the xyz-file {{{

while [ $sz_fo -gt 0 ]; do 
		
		  	sz_fo=` wc -l $fo | awk '{ print $1 }' ` 
			start=$((48*$i_frame-45))
			end=$((48*$i_frame))
			
			gawk 'NR==start, NR==end { print $2, $3, $4 } ' start=$start end=$end $f >& $xyz

			cat $fo | sed -n '1,48p' >& $this_frame
			cat $fo | sed '1,48d' >& $fo.1; mv $fo.1 $fo

			xyz_header=` gawk 'NR==2 { print $0 }' $this_frame ` 

			force=` echo "$xyz_header" | awk ' { print $4 } ' | sed 's/;$//' ` 

			x=( ` cat $xyz | awk '{ print $1 }' ` )
			y=( ` cat $xyz | awk '{ print $2 }' ` )
			z=( ` cat $xyz | awk '{ print $3 }' ` )

			xm=0.0; ym=0.0; zm=0.0

			for i in {0..45}
				do
			  		xm=$( echo "$xm+${x[$i]}" | bc )
			  		ym=$( echo "$ym+${y[$i]}" | bc )
			  		zm=$( echo "$zm+${z[$i]}" | bc )
			done
			  
			xm=$( echo "$xm/$natoms" | bc )
			ym=$( echo "$ym/$natoms" | bc )
			zm=$( echo "$zm/$natoms" | bc )

			srgyr=0.0

		for i in {0..45}
			do
			  dx=$( echo "${x[$i]}-$xm" | bc )
			  dy=$( echo "${y[$i]}-$ym" | bc )
			  dz=$( echo "${z[$i]}-$zm" | bc )
			  dr2=$( echo "$dx^2+$dy^2+$dz^2" | bc )
			  dr2=$( echo "$dr2/$natoms" | bc )
			  srgyr=$( echo "$srgyr + $dr2" | bc ) 
		done
		
			rgyr=$( echo "scale=10; sqrt($srgyr)" | bc )

			eh "File: $f Frame: $i_frame RGYR: $rgyr"
			eh "XYZ comment line: $xyz_header"
			eh 
			eh "Output rgyr file: $rgyr_file"

			echo "$i_frame $force $rgyr" >> $rgyr_file

			i_frame=$(($i_frame+1))
done

		eh "For $1.xyz, number of frames: $i_frame"
		eh "Finished writing output rgyr file: ` basename $rgyr_file` "
		# }}}

