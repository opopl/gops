
# Extract energies from the file lowest, and place them into the output file ef.dat

low_en=( ` grep 'Energy' $dn/lowest | sed 's/^[a-zA-Z 0-9]*=//; s/f[a-z 0-9]*$//' | head -$n_en ` )

gm_energy="${low_en[0]}"

eo "GM energy is = $gm_energy "

#if [ [ $gmin_converged -eq 1 ] || [ $gmin_conv_always -eq 1 ] ]; then 

echo "$var_string ${low_en[@]}" >> $output_dir/ef.dat
