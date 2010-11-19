
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

# copy files to common directories

outd=` basename $output_dir `

eo "Copying: ef.dat to ef/ef_$outd.dat"
cp $output_dir/ef.dat ef/ef_$outd.dat

eo "Appending: ef.dat to ef/ef_all.dat"
cat $output_dir/ef.dat >> ef/ef_all.dat

eo "Copying: gm.xyz to xyz/$outd.xyz"
cp $output_dir/gm.xyz xyz/$outd.xyz

eo "Copying: r.log to rlog/$outd.rlog"
cp $output_dir/r.log rlog/$outd.rlog

eo "scp to $leonov..."

remote_dir="$leonov:~/p46/data/"

scp ef/ef_$outd.dat $remote_dir/ef/
scp xyz/$outd.xyz $remote_dir/xyz/
scp rlog/$outd.rlog $remote_dir/rlog/
