#/bin/bash

#func=`echo $2 | tr '[a-z]' '[A-Z]'`

#files=( geopt.f ) 
funcs=( "GEOM_TO_VAR" "VAR_TO_GEOM" "CHAINBUILD" "UNRESINIT" "INT_FROM_CART" \
"ZEROGRAD" "ETOTAL" "GRADIENT" "LUDCMP" "LUBKSB" ) 
#files=( intbfgsts.f intsecdiag.f )
#files=( keyword.f mylbfgs.f )
#files=( intbfgsts.f )
#files=( path.f unrescalcdihe.f unresconnectsections.f ) 
#files=( ` find . \( -name "*.f" -o -name "*.f90" \)` )
#files=( OPTIM.F )
files=( diis.f intbfgsts.f unresoptim.f )

for file in ${files[@]}; do
for func in ${funcs[@]}; do
	echo $file $func
	cat $file | sed "s/^[ \t]*\(CALL\)[ \t]*\($func\)/!\1 \2/" > n; mv n $file
done
done
