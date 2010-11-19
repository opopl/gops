
#n1=$(($n%2)) 
#n2=$(($n/2)) 

#case "$n2" in
	#0) fr=( 0.0001 0.001 ) ;;
	#1) fr=( 0.001 0.01 ) ;;
	#2) fr=( 0.01 0.1 ) ;;
	#3) fr=( 0.1 1.0 ) ;;
#esac

#case "$n1" in
          #0) sys="G46" ; od="go$n2" ; queue=s1 ;;
	#1) sys="P46" ; od="wt$n2" ; queue=s1 ;;
#esac

#let nf=9

radius=2.0; nrg=1 ; temp=0.03 ; nst=2; queue=s2

fr=( 0.0 0.0001 )

let nf=100

sys="P46"
nsteps=100000

force_min="${fr[0]}"
force_max="${fr[1]}"
force_i=$( echo "scale=10; ($force_max-$force_min)/$nf" | bc )

od="$od"_N_`printf "%1.0e\n" $nsteps`
od="$od"_f_`printf "%1.1e\n" $force_min`
od="$od"_`printf "%1.1e\n" $force_max`

od="out"

source $shd/define_forces_0.sh 

#case $nst in
          #3) nsteps=1000000 ;;
          #2) nsteps=100000 ;;
          #1) nsteps=10000 ;;
          #0) nsteps=1000 ;;
#esac

