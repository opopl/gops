
radius=2.0; nrg=1 ; temp=0.03 ; nst=5; queue=s2
sys="G46"

case "$sb_mode" in
  	1) # {{{
		n1=$(($n%2)) 
		n2=$(($n/2)) 
		
		case "$n2" in
			0) fr=( 0.0001 0.001 ) ;;
			1) fr=( 0.001 0.01 ) ;;
			2) fr=( 0.01 0.1 ) ;;
			3) fr=( 0.1 1.0 ) ;;
		esac
		
		case "$n1" in
			0) sys="G46" ; od="go$n2" ; queue=s1 ;;
			1) sys="P46" ; od="wt$n2" ; queue=s1 ;;
		esac
		
		let nf=9
		;;
	      # }}}	
	2) fr=( 0.0 0.0001 ) ; nst=4; let nf=100 ; od="out" ;;
esac

force_min="${fr[0]}"
force_max="${fr[1]}"
force_i=$( echo "scale=10; ($force_max-$force_min)/$nf" | bc )
nsteps=$( echo "10^$nst" | bc )

#od="$od"_N_`printf "%1.0e\n" $nsteps`
#od="$od"_f_`printf "%1.1e\n" $force_min`
#od="$od"_`printf "%1.1e\n" $force_max`

source $shd/define_forces_0.sh 


