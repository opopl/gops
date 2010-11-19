
eo "Run mode: $mode"
eo
eo "Output data directory: $output_base_dir"
eo

case "$mode" in
  	"-all")
		source $shd/define_forces.sh
		source $shd/adjust_data.sh
		
		eo "$prog input file, after sourcing adjust_data.sh : "
		eo "=============================="
		cat data.G46 >> $main_log_file
		eo "=============================="
		
		source $shd/reset_files.sh
		source $shd/eo_pars.sh
		
		source $shd/loop_forces.sh
		source $shd/end.sh
	;;
esac

