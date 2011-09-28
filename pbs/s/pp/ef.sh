#/bin/bash

# directory where this script resides
export shd="`dirname $(readlink -f $0)`"

# name of this script 
export this_script=` basename $0 `

deps=( ef ea ea_fz ea_fvsz ea_ta ea_ba ea_nbond ea_bond ) 
rm -rf *.ps *.pdf
all="a.ps"

while [ ! -z "$1" ]; do 
          #intro {{{
	for dep in ${deps[@]}; do
		xlabel="f"

		sys="$1" 
		depname=`echo $dep | sed 's/_.*$//g'`
		main_dfile="$depname"."$sys".sh
		dfile=$main_dfile
		plfile="$depname.$sys.ps"

		case $sys in
	  		G46) stitle="Go-like" ;;
	  		P46) stitle="Wild-type" ;;
		esac
	  #}}}
		#dep switch {{{
		xlabel="Force f"
		case $dep in 
	  		ef) xn=1; yn=2; 
					btitle="Total Energy vs Force" 
					ylabel="E"
					;;
	  		ea_fz) xn=1; yn=7; 
					btitle="-F*dZ vs Force" 
					ylabel="-F*dZ"
					;;
			ea_fvsz) 
					xn=1; yn=2
					dfile="ea_fvsz.sh"
					dfz.pl $main_dfile > $dfile
					btitle="Extension vs Force" 
					xlabel="Force f"
					ylabel="Extension |dZ|"
					;;
			ea_ta) xn=1; yn=6; 
					btitle="Torsional Angle Energy vs Force" 
					ylabel="E"
					;;
			ea_ba) xn=1; yn=5; 
					btitle="Bond Angle Energy vs Force" 
					ylabel="E"
					;;
			ea_nbond) xn=1; yn=3; 
					btitle="Non-Bonded Energy vs Force" 
					ylabel="E"
					;;
			ea_bond) xn=1; yn=4; 
					btitle="Bonded Energy vs Force" 
					ylabel="E"
					;;
		esac
	    #}}}
	#gnuplot {{{
gnuplot << gp
	set term postscript enhanced
	set output "$plfile"	
	set title "$btitle, $stitle"
	set xlabel "$xlabel"
	set ylabel "$ylabel"
	plot "$dfile" using $xn:$yn with points
gp
#}}}

		cat $plfile >> $all
		echo "false 0 startjob pop" >> $all
	
	done 
	shift
done
		
ps2pdf $all
