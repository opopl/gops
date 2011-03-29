#/bin/bash

#func=`echo $2 | tr '[a-z]' '[A-Z]'`

files=( geopt.f ) 
funcs=( "GEOM_TO_VAR" "VAR_TO_GEOM" "CHAINBUILD" ) 

for file in ${files[@]}; do
for func in ${funcs[@]}; do
	echo $file $func
	cat $file | sed "s/^[ \t]*\(CALL\)[ \t]*\($func\)/!\1 \2/" > n; mv n $file
done
done
