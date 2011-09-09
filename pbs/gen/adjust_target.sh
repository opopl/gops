
if [ $gmin_conv_always -eq 1 ]; then 
  	case "$adj_target" in
	  	"1") target_energy="$gm_energy"
	       	     eo "Target energy reset to the current GM energy"
		     ;;
	esac
fi
