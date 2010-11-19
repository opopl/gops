
let i=0
for force in "${forces[@]}"
	do
	  	echo "$force"
		$shd/r -q $queue -nrg $nrg -od $od -radius $radius -nsteps $nsteps -temp $temp -ap -force $force 
		sleep $n_sleep
		i=$(($i+1))
done
