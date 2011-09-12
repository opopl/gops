
echo "	RUNS:"
	let i=0
	for rn in "${runs[@]}"
		do
		  echo "$i	$rn"
		  i=$(($i+1))
	done

