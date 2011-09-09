
F=$force_min
let i=0
let fabove=0

while [ "$fabove" -eq "0" ]
	do
		forces[$i]=$F;
		F=$( echo "$F + $force_i" | bc );
		fabove=$( echo "$F > $force_max" | bc )
		
		let i=$(($i+1));
done

