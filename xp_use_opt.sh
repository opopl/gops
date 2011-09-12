
#ln=( ` cat $script_dir/opt.txt | sed -n "/$mod/p" ` )


#echo "Using opt.txt..."
#echo "	Values for $mod: "
#echo "		xamin $xamin"
#echo "		yamin $yamin"
#echo "		zamin $zamin"

case "$file0" in
	"chf") 
		case $i_frame in
		  	1 | 3 | [4-9] | 1[0-9]) 
			ln=( 45 90 45 100 ) 
			;;
			2) 
			ln=( 90 90 45 70 )
			;;
			20) 
			ln=( 60 0 70 70 ) 
			;;
			21) 
			ln=( 50 54 0 70 )
			;;
		esac
	;;
esac	

xamin=${ln[0]} ; yamin=${ln[1]} ; zamin=${ln[2]}
xamax=$xamin ; yamax=$yamin ; zamax=$zamin
zoom=${ln[3]}
