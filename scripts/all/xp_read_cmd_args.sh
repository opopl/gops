#!/bin/bash - 
#===============================================================================
#
#          FILE:  xp_read_cmd_args.sh
# 
#         USAGE:  ./xp_read_cmd_args.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 10/11/10 13:04:13 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

while [ ! -z "$1" ]; do
	  	case "$1" in
		  	# ==================
			#
		  	# Ground state of the WT/Go model
			#
		  	"-0") 
				ax=$2; shift; ay=$2; shift; az=$2; shift; shift
				./xp -ax $ax -ay $ay -az $az \
					--in_file gogs \
					--in_dir ../x0 \
					--pics_dir ../pics-gogs \
				       	--picoutfile gogs-$ax-$ay-$az $* ; 
				exit 
				;; 
			"-00") ./xp -0 45 90 45 $* ; exit ;;
			"-0x") shift; ./xp -00 -pof gogs-rx -rx 0 180 $*; exit  ;;
			"-0y") shift; ./xp -00 -pof gogs-ry -ry 0 180 $*; exit  ;;
			"-0z") shift; ./xp -00 -pof gogs-rz -rz 0 180 $*; exit  ;;
			#
			# printing text in output images
			#
			"-nt") print_text="no" ;;
			"-pra") print_text="ra" ;;
			#
			# ===================
			#
		  	"-rm") rm -f $pics_dir/*.$pic_ext ; exit ;;
	    		"-f" | "--in_file" ) file0="$2" ;;
			"-p" | "--pics_dir") pics_dir="$2"
				mkdir -p $pics_dir ;;
			"-zm") zoom=$2 ;;
			"-opt") use_opt=1 ;;
			"-rot") rotate=1; no_rotate=0 ;;
			"-pof" | "--picoutfile") reset_pof=1 ; pof="$2"."$pic_ext" ;;
			"-a") source "$this_script"_cmd"$1".sh ;;
			"-aa") append_angles=1 ;;
			"-af") append_force=1 ;;
			"-naf") append_force=0 ;;
			"-ax") 
				xamin=$2 
				xamax=$xamin
			;;
			"-ay") 
				yamin=$2 
				yamax=$yamin
			;;
			"-az") 
				zamin=$2
				zamax=$zamin
			       	;;

			"-rx") rx=1; no_rotate=0; xamin=$2; xamax=$3; xai=$4 ;;
			"-ry") ry=1; no_rotate=0; yamin=$2; yamax=$3; yai=$4 ;;
			"-rz") rz=1; no_rotate=0; zamin=$2; zamax=$3; zai=$4 ;;

			"-d" | "--in_dir") xyz_dir="$2" ;;
			"-l") ls $xyz_dir ; exit ;;
			"-fmin") frame_min=$2 ; reset_f=1 ;;
			"-fmax") frame_max=$2 ; reset_f=1 ;;
			"-fr") frame_min=$2 ; frame_max=$2 ; reset_f=1 ;;
			"-fns") font_size=$2 ;;
			"-fa") reset_f="2" ;;
			"-log") write_log_file=1 ;;
	  	esac
	  	shift
	  done

