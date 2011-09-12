
# Extract energies from the file lowest, and place them into the output file ef.dat

low_en=( ` grep 'Energy' $dn/lowest | sed 's/^[a-zA-Z 0-9]*=//; s/f[a-z 0-9]*$//' | head -$n_en ` )

gm_energy="${low_en[0]}"
gme=$( echo "scale=10; ($gm_energy-($gs_energy))/($gs_energy)" | bc )

eo "GM energy is = $gm_energy "
eo "Number of extracted lowest energies: ${#low_en[@]}"

outf=$output_dir/ef.dat

case "$print_ef" in
  	"a") echo "` printvar 'force_long' $force` $gme ${low_en[@]}" >>  $outf ;;
#echo "$var_string ${low_en[@]}" >> $output_dir/ef.dat
esac
