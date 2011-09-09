
end_date=`  date_hms `
end_time=`  date_in_secs `

pref="end"

eo "end date=$end_date "
eo "end time=$end_time "


# remove blank lines in output files

for of in ${output_files[@]}
	do
	  sed '/^$/d' $of > n
	  mv n $of
done

total_time=` time_dhms $end_time $start_time `

eo "total time elapsed= $total_time " 

echo "TIME: $total_time"

