#!/bin/bash 

export shd="`dirname $(readlink -f $0)`"
source $shd/f.sh

#for n in {0..7}
	#do
		  #source sub_sh_set_pars.sh
		#opts0="-q $queue -sys $sys -nrg $nrg -od $od -radius $radius -nsteps $nsteps -temp $temp -nap"
		#opts="$opts0 -if $force_min $force_max $force_i"
		#r $opts
		#sleep 3
#done
# dirs: go*, wt*

source $shd/sub_sh_set_pars.sh

opts0="-q $queue -sys $sys -nrg $nrg -od $od -radius $radius -nsteps $nsteps -temp $temp -nap"

opts="$opts0 -if $force_min $force_max $force_i"

$shd/r $opts
