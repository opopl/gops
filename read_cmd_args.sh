
export logf=rlog.tmp
touch $logf

eo "reading commandline arguments..." 

cmd_args=( ` cat cmd_args ` )

n_cmd_args=${#cmd_args[@]}

if [ $n_cmd_args -eq 0 ]; then
	eo "... no commandline arguments supplied"
else
	eo "Number of command-line arguments: $n_cmd_args"
	eo "Command-line arguments: ${cmd_args[@]}"
	let i=0
for arg in ${cmd_args[@]} # {{{
	do
	  arg2=${cmd_args[$(($i+1))]}
	  arg3=${cmd_args[$(($i+2))]}
	  arg4=${cmd_args[$(($i+3))]}
		case "$arg" in
		  	"-ap" | "--append_files")  append_gm=1 
				append_files=1
				;;
			"--append_job-name") append_job_name=1 ;;
			"-nap" | "--not_append_files") append_files=0  ;;
			"--plot_ef" | "-pef" | "--plot_gm_xyz" | "-pgm" )
		       		mode="$arg"	
				append_files=1 
				;;
			"-pstart") 	pull_start=$arg2 ;;
			"-pend") 	pull_end=$arg2 ;;
			"-nsteps") 	nsteps=$arg2 ;;
			"-te") 		target_energy=$arg2 ;;
			"-temp") 	temp=$arg2 ;;
			"-n_en") 	n_en=$arg2 ;;
			"-force") 	force=$arg2 ; get_force="single" ;;
			"-nrg") 	nrg=$arg2 ;;
			"-radius") 	radius=$arg2 ;;
			"-seed")	seed=$arg2 ;;
			"--no_force_loop") do_force_loop=0 ;;
			"--force_loop") do_force_loop=1 ;;

			# PBS, for information purposes
			
			"-q" | "--pbs_queue_name") queue_name="$arg2" ;;
			
			# output directory

			"-od" | "--output_dir_name") 
				export output_base_dir="$arg2" 
				export output_dir="$data_dir/$output_base_dir" 
			  	mkdir -p $output_dir
				;;

			# JOBS

			"--same-output-files")  ;;
			"-j") job_number=$arg2 ;;
			"--append_job_number") append_job_number=1 ;;
			"--n_append_job_number") append_job_number=0 ;;

			# GMIN

			"--print_gmin_warnings") print_gmin_warnings=1 ;;			
			"--print_acceptance_ratios") print_acceptance_ratios=1 ;;
			"--remove_dn") remove_dn=1 ;;
			"--not_remove_dn") remove_dn=1 ;;

			"-bs" | "--basin-sampling") ;;
			
			"-sys") sys="$arg2" ;;

			# OPTIM 

			# ...> 

			"--pref_pbs_job_id") pref_pbs_job_id=1 ;; 
			"--pref_job_name") pref_job_name=1 ;; 
		
			# PULLING
	
			"--adj_target") adj_target="$2" ;;
			"-if") 
				get_force="increment" 
			      	do_force_loop=1

				force_min=$arg2
				force_i=$arg4
				force_max=$arg3
				;;
			"-af")  get_force="array" ;;
		esac
		let i=$(($i+1))
done # }}}
fi

export main_log_file=$output_dir/r.log

#if [ -f "pre.log" ]; then # {{{
	#if [  $append_files -eq 1 ]; then 
			#cat "pre.log" >> $main_log_file
		 #else
			#cat "pre.log" >& $main_log_file
	#fi	  
#fi # }}}

cat $logf >> $main_log_file

rm $logf
rm -f pre.log

export logf=$main_log_file
