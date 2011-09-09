
F=$force_min
let i=0
let fabove=0

while [ "$fabove" -eq "0" ]
do
forces[$i]=$F;
F=$( echo "$F + $force_i" | bc );
fabove=$( echo "$F > $force_max" | bc )
#eo "fabove=$fabove; F=$F; force_min=0.0; force_max=$force_max; force_i=$force_i"

let i=$(($i+1));
done

