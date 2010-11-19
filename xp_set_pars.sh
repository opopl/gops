#!/bin/bash - 
#===============================================================================
#
#          FILE:  xp_set_pars.sh
# 
#         USAGE:  ./xp_set_pars.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 10/11/10 13:10:16 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

r="$shd"

xyz_dir="$r/../xyz/"
tmp_dir="$r/tmp"
pics_dir="$r/../pics/"
file0="chf"

reset_f=0
reset_pof=0

frame_min=1
frame_max=1
font_size=30
show_date=0 
write_log_file=1
pic_ext="gif"
rotate=0
no_rotate=1
use_opt=0

append_angles=0
append_force=0
print_text="all"

zoom=70

rx=0 ; ry=0 ; rz=0

xai=10 ; yai=10 ; zai=10
xamin=0; yamin=0; zamin=0
xamin=45; yamin=45; zamin=45

