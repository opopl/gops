#!/bin/bash - 
#===============================================================================
#
#          FILE:  xp_print_help.sh
# 
#         USAGE:  ./xp_print_help.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 10/11/10 13:33:16 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

cat << EOF
================================================================================
SCRIPT NAME: $this_script
PURPOSE:
	convert a given XYZ file to an image (gif) file
USAGE:
	$this_script OPTS
OPTS:
	============

	SPECIAL

	-0 AX AY AZ OPTS 	=> xp -ax AX -ay AY -az AZ 
					--in_file gogs
					--in_dir ../x0 
					--pics_dir ../pics-gogs 
					-pof gogs-ax-ay-az OPTS
				
					Input XYZ file: gogs.xyz
					Input XYZ directory: ../x0
					Output directory for images: ../pics-gogs
					Output image filename: ax-ay-az
						with ax, ay, az the rotation angles.
	-0x INC 		=> xp -0 -rx 0 180 INC 			
	-0y INC 		=> xp -0 -ry 0 180 INC 			 
	-0z INC 		=> xp -0 -rz 0 180 INC 			 

	============
	TEXT IN IMAGES

	-nt				no text in images
	-pra				print only rotation angles 
	
	COMMON
	

	-rm 				= rm -f $pics_dir/*.$pic_ext
	-f, --in_file file0			xyz filename (without xyz-extension!)
	-a				all xyz files in a given xyz-directory, all frames
	-d,--in_dir xyz_dir			specify directory for xyz files
	-p, --pics_dir pics_dir			directory for images
	-pof, --picoutfile NAME		pof name for the output image file (specify without extension!)
	-l 				list xyz files

	============
	FRAMES
	
	-fmax		frame_max		maximal frame number
	-fmin		frame_min		minimal frame number
	-fr		frame_index		only this frame
	-fa		all frames in the given file(s)
	
	============
	
	-fns		font_size		font size
	-log		write output into log-file
	-rot		rotate x, y, z angles

	------------
	-opt		use 'opt.txt' file which contains data for 
				optimal parameters.
				Format for file opt.txt:
					FILE_LABEL xa ya za zoom   

			# means comment
	------------

	-rx 		rotate x angle
	-ry		rotate y angle
	-rz		rotate z angle

	-ax		x angle
	-ay		y angle
	-az		z angle
	-aa		append angles info to the output image files
	-af		append forces to the output image file
	-naf		overwrite image file for the same value of force
	-zm	zoom 	zoom value, default is $zoom

DEFAULTS:
	xyz_dir=$xyz_dir
	pics_dir=$pics_dir
	file0=$file0
	frame_min=$frame_min
	frame_max=$frame_max
	font_size=$font_size
	write_log_file=$write_log_file
	pic_ext="$pic_ext"
	append_force=$append_force

	------------------
	Rotations:
	------------------
	rotate=$rotate
	rx=$rx,	ry=$ry,	rz=$rz
	xai=$xai, yai=$yai, zai=$zai
	xamin=$xamin, yamin=$yamin, zamin=$zamin
	xamax=$xamax, yamax=$yamax, zamax=$zamax
	use_opt=$use_opt
	append_angles=$append_angles
	zoom=$zoom
================================================================================
EOF

