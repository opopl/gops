
eo "entering force iteration loop..." 

pref_old="$pref"

export pref=$pref."force_iterate"

for force in "${forces[@]}"
	do
		dn="$data_dir/dnj/dn-j$pbs_job_id"

		eo "force is now = ` printvar 'force' $force ` "
		eo "test directory is = $dn "

		echo "$dn" >& $data_dir/current_dn

		$shd/adjust_data.sh 
		eo "Running $exe..."
		
		start_time_1=` date_in_secs `

		$exe $radius $dn $seed $nrg  

		cp $dn/lowest xyz/lowest/f-$force-j$pbs_job_id.xyz
		
		var_string="$force $temp $radius $nrg $nsteps"
		#var_full_string="` date_dm_hm` f=`printvar 'force' $force` T=$temp R=$radius nrg=$nrg N=$nsteps j=$pbs_job_id s=$sys"
		var_full_string="` date_dm_hm` f=`printvar 'force' $force`"

		source $shd/extract_energies_from_lowest.sh
		source $shd/adjust_target.sh
		source $shd/append_gm.sh
		source $shd/append_gmin_out.sh
		source $shd/eo_warnings.sh
		source $shd/append_hits.sh
		source $shd/append_pdf.sh
		source $shd/append_energy_markov_dat.sh

		if [ $remove_dn -eq 1 ]; then 
			rm -rf $dn
	      	fi
		
		end_time_1=` date_in_secs `
		calc_time_1=` time_dhms $end_time_1 $start_time_1 `
		calc_time_1_secs=$(($end_time_1-$start_time_1))

		eo "calculation time: $calc_time_1, $var_full_string"

		if [ $write_t_dat -eq 1 ]; then  
			timing_data=" ` date_dm_hm` $var_string $calc_time_1_secs"
			echo "$timing_data $gm_energy" >> $data_dir/t.dat
			echo "$timing_data $gm_energy" >> $output_dir/tl.dat
	        fi

		# }}}
done

export pref="$pref_old"
eo "finished force iterate loop"
