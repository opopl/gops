
rm -f $a.pdf
let i=0

for data_file in "${data_files[@]}"
  	do
	  fs=` basename $data_file `
	  # now, for fs, get rid of the prefix ef_ or similar one ...
	  fs=` echo $fs | sed "s/^$prefix//" ` 
	  # finally, get rid of the dat extension
	  fs=${fs%.dat}

	  #title=" ` set_title $fs `"
	  file_header=` awk 'NR==1' $data_file  `
	  # 

	  input_file_date=` echo $file_header | awk -F ";" '{ print $1 }' | sed 's/^[#@@]*//'` 
	  temp=` echo $file_header | awk -F ";" '{ print $4 } ' | sed 's/^[ T=]*//'` 
	  nsteps0=` echo $file_header | awk -F ";" '{ print $2 } ' | sed 's/^[ N=]*//'` 
	  radius=` echo $file_header | awk -F ";" '{ print $6 } ' | sed 's/^[ R=]*//'` 
	  sys=` echo $file_header | awk -F ";" '{ print $8 }' | sed 's/^[ sys=]*//'` 
	  
	  
	  output_file="$a$i.ps"
	  title=""

	  title="$title Input file date: $input_file_date \n"

	  if [ $print_date -eq 1 ]; then
		    	title="$title Plotting Time: `date_dm_hm` \n"
	  fi

	  fs=` echo $fs | sed 's/_//'`

	  if [ $print_file -eq 1 ]; then
		    title="$title Input File: $fs \n"
	  fi
	  
	  title="$title T=$temp, N=$nsteps0,  R=$radius; System $sys \n" 
	 
	  source ../gp_cat_gnuplot.sh

echo "Data type: $prefix"
echo "Data file: $data_file, plot file: $output_file, its size: ` wc -c $output_file | awk '{ print $1 }' ` "
echo "Plot title: $title"

ps2pdf $output_file
pdf_add $a.pdf $a$i.pdf
i=$(($i+1))
#
#
################################
done

echo "Final file: $a.pdf, its size: `wc -c $a.pdf | awk '{ print $1 }'`"

cp $a.pdf $plot_dir
scp $a.pdf $venus:~/
scp $a.pdf $leonov:~/p46
rm *.ps *.pdf

