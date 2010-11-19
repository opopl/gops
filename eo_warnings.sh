
warning_words=( "Quench" )

for warning_word in "${warning_words[@]}"
	do
		failed_quench_numbers= \
		( ` cat $dn/GMIN_out | grep 'WARNING' |  \
			awk '/ww/' ww="$warning_word" | \
			awk -F "-" '{ print $2 }' |  \
			sed 's/[a-z ]*$//' ` ) 
done

number_of_failed_quenches=${#failed_quench_numbers[@]}
nfq=$number_of_failed_quenches

if [ $nfq -gt 0 ]; then
	 eo "Warnings:   "
	 echo "==========="
 	 eo "Number of failed quenches: $nfq"	
	 eo "Failed quenches are: ${failed_quench_numbers[@]}"
fi 

cat $dn/GMIN_out | grep 'Final Quench' | tail -2  >> $logf
echo "==========================" >> $logf
