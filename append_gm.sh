
# Append the GM conformation to a separate file called gm.xyz

tmpxyz="tmp/gm-temp.xyz"
head -1 $dn/lowest > $tmpxyz
#
# 	In the temporary XYZ file, which is called here $tmpxyz, we will need to include information about: 
#		1) global minimum energy 
#		2) value of the force
# 		3) maybe other parameters (like time/date etc.?) if necessary
#

gm_energy_s="GM energy= ${low_en[0]}; "
gm_energy_s="E= ${low_en[0]}; "
force_s="force= $force; "
force_s="F= $force; "
steps_s="MC steps= $nsteps; "
steps_s="N= $nsteps; "

other_s="  ` date_dmy_hms ` ; nrg=$nrg; T=$temp; R=$radius" 
xyz_comment_string="$gm_energy_s $force_s $steps_s $other_s"

# Replace the comment (the second one) line in $tmpxyz  with a new one:
echo "$xyz_comment_string" >> $tmpxyz
sed -n '3,48 p' $dn/lowest >> $tmpxyz

# Now we may append the tmpxyz file to the gm.xyz
cat $tmpxyz >> $output_dir/gm.xyz
