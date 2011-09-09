
eo "entering force iteration loop..." 

pref_old="$pref"

export pref=$pref."force_iterate"

for force in ${forces[@]}
	do
		#dn=f-$f
		mkdir -p dnj
		dn="dn-j$pbs_job_id"

		eo "force is now = $force"
		eo "test directory is = $dn "

		echo "$dn" >& current_dn

		change.sh force $force
		change.sh target $target_energy

		eo "Target energy: $target_energy "
		eo "Running $exe..."
		
		start_time_1=` date_in_secs `

		$exe $radius $dn $seed $nrg  

		cp $dn/lowest xyz/lowest-f-$force.xyz
		
		if [ $print_gmin_warnings -eq 1 ]; then 	
			grep $dn/GMIN_out 'WARNING' >> $logf  
	 	fi
		
		var_string="$force $temp $radius $nrg $nsteps"
		var_full_string="f=$force T=$temp R=$radius nrg=$nrg N=$nsteps j=$pbs_job_id"

		source extract_energies_from_lowest.sh
		source adjust_target.sh
		source append_gm.sh
		source append_gmin_out.sh
		source append_hits.sh
		source append_pdf.sh
		source append_energy_markov_dat.sh

		#source append_acceptance_ratios.sh
		if [ $remove_dn -eq 1 ]; then 
			rm -rf $dn
	      	fi
		
		end_time_1=` date_in_secs `
		calc_time_1=` time_dhms $end_time_1 $start_time_1 `
		calc_time_1_secs=$(($end_time_1-$start_time_1))

		eo "calculation time: $calc_time_1, $var_full_string"

		if [ $write_t_dat -eq 1 ]; then  
			timing_data=" ` date_dm_hm` $var_string $calc_time_1_secs"
			echo "$timing_data $gm_energy" >> t.dat
			echo "$timing_data $gm_energy" >> $output_dir/tl.dat
	        fi

		# }}}
done

export pref="$pref_old"
eo "finished force iterate loop"
