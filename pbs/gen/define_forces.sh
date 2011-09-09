
eo "define_forces.sh: define the array of forces"

case "$get_force" in 
		"increment")
			eo "get_force=$get_force" 
			eo "force_min=$force_min, force_i=$force_i, force_max=$force_max" 
			eo
			
			F=$force_min
			let i=0
			fabove=0

			while [ "$fabove" -eq "0" ]
			   		do	
						forces[$i]=$F;
						F=$( echo "$F + $force_i" | bc );
						fabove=$( echo "$F > $force_max" | bc )
						#eo "fabove=$fabove; F=$F; force_min=0.0; force_max=$force_max; force_i=$force_i"

						let i=$(($i+1));
			done
		;;
		"array") forces=( 0.0 0.01 0.02 0.1 0.5 1.0 ) ;;
		"single") forces=( $force ) ;;
esac

if [ $do_force_loop -eq 0 ]; then 
		forces=$force
fi	

# number of sampled force points

nforces=` echo ${#forces[@]} `
#eo "$nforces"
