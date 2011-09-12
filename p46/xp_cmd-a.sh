#!/bin/bash - 
#===============================================================================
#
#          FILE:  xp_cmd-a.sh
# 
#         USAGE:  ./xp_cmd-a.sh 
# 
#   DESCRIPTION:  G
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 10/11/10 13:18:33 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

echo "XYZ directory: $xyz_dir"

cd $xyz_dir

for xyz_file in ` ls *.xyz `
	do
		xyz_bare_file=${xyz_file%.xyz}	
		echo "Processing $xyz_file... "

		sz=$(stat -c%s "$xyz_file")
		echo "File size: $sz"

		if [ $sz -gt 0 ]; then
			$this_script -f "$xyz_bare_file" -fa
		fi
done
cd -

