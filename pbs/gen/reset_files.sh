
let i=0

for f in ${output_files[@]}
		do
		  	if [ $append_files -eq 0 ]; then
					echo ${output_file_lines[$i]} >& $f
			      else
					echo ${output_file_lines[$i]} >> $f
			fi
			let i=$(($i+1))
done
