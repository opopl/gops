#!/bin/bash

source f.sh

this_script_dir="` dirname $(readlink -f $0)`"
this_script=` basename $0 `
tsd=$this_script_dir
plot_dir=$tsd/plots

out_dirs=( ef rgyr )
out_dirs=( ef )
ef_dir="ef"

print_file=1
print_date=1

if [ -z $1 ]; then
cat << EOF 
====================================================
NAME: gp
PURPOSE: plot data files in ${out_dirs[@]}
OPTS:
	a		all in ${out_dirs[@]}
	-ef NAME 	plot only file ef_NAME.dat		
	-pd		print date
====================================================
EOF
	else
	  while [ ! -z $1 ]; do 
		case "$1" in
		 	a)  source gp_all.sh ;;
			-ef) 
				a=$2 
				prefix="ef"

			cd $ef_dir
				data_files=( "$prefix"_$a.dat )
				source $this_script_dir/gp_loop_over_dat.sh
			cd $this_script_dir

			;;
			-pd)
			print_date=1
			;;
		esac
		shift
	done
################################
fi
