
eo "Run mode: $mode"
eo
eo "Output directory: $output_dir"
eo

case "$mode" in
  	"-all")
source define_forces.sh
source adjust_data.sh
source reset_files.sh
source eo_pars.sh
source loop_forces.sh
#source plot_gm_xyz.sh
#source plot_e_vs_force.sh
source end.sh
	;;
	"-plot_ef" | "-pef") 		source plot_e_vs_force.sh ;;
	"-plot_gm_xyz" | "-pgm")	source plot_gm_xyz.sh ;;
esac

