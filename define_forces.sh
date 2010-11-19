
eo "define_forces.sh: define the array of forces"

case "$get_force" in 
		"increment")
			eo "get_force=$get_force" 
			eo "force_min=` printvar 'force' $force_min`" 
			eo "force_max=` printvar 'force' $force_max`" 
			eo "force_i=` printvar 'force' $force_i`" 
			eo
			
			F=$force_min
			let i=0
			fabove=0

			while [ "$fabove" -eq "0" ]
			   		do	
						forces[$i]=$F;
						F=$( echo "$F + $force_i" | bc );
						fabove=$( echo "$F > $force_max" | bc )
						let i=$(($i+1));
			done
		;;
		"single") forces=( $force ) ;;
esac

if [ $do_force_loop -eq 0 ]; then 
		forces=$force
fi	

nforces=` echo ${#forces[@]} `
