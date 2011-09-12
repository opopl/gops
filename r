#!/bin/bash

export shd="`dirname $(readlink -f $0)`"
export shd="/home/op226/p46/sh"
source $shd/sh_header.sh

rm -f cmd_args

cp $data_dir/data.old $shd/data.G46

source $shd/set_base_pars.sh
source $shd/set_pbs_pars.sh
source $shd/set_runs.sh

add_run_opts=""
r_mode="qsub"

echo "$@" > cmd_args

if [ -z "$*" ]; then 
		source $shd/qsub.sh 
      else
		while [ ! -z "$1" ]; do # process cmd arguments {{{ 
			case "$1" in
			  	"-hr") r_mode="runr" ; run_opts="-h runs" ;;
			  	"-hs") r_mode="runr" ; run_opts="-h scripts" ;;
		  		"-h")   r_mode="help"
					help_switch="$2" 
					source $shd/"$this_script"_print_help.sh ;;
				"-q" | "--pbs_queue_name") 	export queue_name="$2" ;;
				"-wt" | "--pbs_wall_time") 	export wall_time="$2" ;;
				"--pbs_job_name")   		export job_name="$2" ;;
				# RUNS
				"-r")  
					run_index="$2" 
					run_opts="${runs[$run_index]}"	
					shift; shift;
					add_run_opts="$@"
					echo "========================" >& pre.log
					echo "RUN: $run_index" >> pre.log
					echo "RUN CMDS: $run_opts  " >> pre.log
					echo "ADDITIONAL CMDS: $add_run_opts  " >> pre.log
					echo "START DATE: ` date_dm_hm ` " >> pre.log
					r_mode="runr"
			       		;;	
				*)	;;
			esac
		shift
	        done # }}}

		case "$r_mode" in
		 	"qsub") source $shd/qsub.sh ;;
			"runr") $shd/r $run_opts $add_run_opts ;;
		esac
fi

# vim:set ft=sh
