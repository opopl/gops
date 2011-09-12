#!/bin/bash
#vim:set fdm=marker

source f.sh
source set_base_pars.sh
source set_gmin_pars.sh

rpdir="$HOME/wrk/rep/trunk/pics/"

view_plot=1
make_plot=1
plot_type="plot"

reset_pf=1

arg="ef"
ext="eps"
gpi="0"
is=2

gp_xsize=1.0 ; gp_ysize=$gp_xsize
gp_font="Helvetica 30"

lx=0.7 ; ly=0.5; lshift=0.05; lk="graph"

colx=1 ; coly=6

# plot title options

echo_date_plotted=1 
echo_out_dir=1 

if [ -z "$1" ]; then
cat << EOF
==============================================================
SCRIPT NAME: pl
PURPOSE: generate plots from output data files
USAGE: pl OPTS
OPTS:
	-d 		 use defaults

	-ps, -eps, -pdf

	-it string	 inner title specified by string

	-gpi			
		0	- 	usual E vs f 
		logxy	-	log dE/E0 vs log f  
	

	GNUPLOT

	-is starting index in the date file, default is 2
	-lx
	-ly
	-lshift
	-lk word		graph, screen, first, second

	-pf FILE		output plot file			

	------------------

	-p arg		 arg specifies what to plot
	-f ifile	 name for specific input file.
				Default convention is that ifile=arg.
				i.e., pl -p ef will try to plot file ef.dat
				this option should be specified AFTER the -p one
	---------
	it can be:
		   ef	-  energy vs force
			   default value is arg=ef
		   gm  	- 
		   	   plot gm.xyz using pymol 
		   pdf  - 
		   statistics for the mean global minimum first encounter time

		   s_f0	- energy.dat, markov.dat for f=0
	---------
	-ac i	-	 use archive output dir with index i
	-v 	-	 view plot file after plotting
	-od	-	 output directory
==============================================================
EOF
	else
		mkdir -p plots
	  	while [ ! -z $1 ]; do 
			case "$1" in
			  	-0) pl -f $2 -pf $2 -gpi logxy -12 ;;
			  	-00) pl -0 $2 -is $3 ;;
			  	-12) colx=1; coly=2;;
			  	-16) colx=1; coly=6;;
				-pf) pf=$2 ;;
			  	-is) is=$2 ;;
			  	-lx) lx=$w ;;
			  	-ly) ly=$2 ;;
			  	-lk) lk=$2 ;;
			  	-lshift) lshift=$2 ;;
			  	-od) output_dir="$2" ;;
			  	-d) arg="ef" ; make_plot=1 ; view_plot=1 ; ext="ps" ;;
			  	-v) view_plot=1 ;;
				-f) ifile=$2 ;;
				-p) arg=$2 
					ifile=$arg
				case "$arg" in
					"gm") plot_type="xyz"  ;;  	
					"s_f0") output_dir="f-0.0" ;;
				esac
				;;
				-it) inner_title=$2 ;; 
				-ac) 	
					let ao_count=$2
					ao_count=$(($ao_count+1)) 
					awk_string="NR==$ao_count { print \$2 }"
			                export output_dir=` awk "$awk_string" ao.log `	
					;;
				-ps) ext="ps" ;;
				-eps) ext="eps" ;;
				-pdf) ext="pdf" ;;
				-gpi) gpi=$2 ;;
			esac
			shift
		done

		if [ $make_plot -eq 1 ]; then  

			if [ $echo_date_plotted -eq 1 ]; then 
				plot_title="Plotted on ` date_dmy_hms `	$plot_title"
		        fi
			if [ $echo_out_dir -eq 1 ]; then 
				plot_title="Data directory: $output_dir  $plot_title"
		        fi

			source plot_arg.sh
		fi

		if [ $view_plot -eq 1 ]; then 
			$psv $gv_opts $plot_file
			#$pdfv $latex_plot_file
		fi
fi
