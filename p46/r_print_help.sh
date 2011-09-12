
cat << EOF
======================================================================
EOF

case "$help_switch" in
  	"" | "all") help_sections=( header exec gmin scripts pbs pars runs opts ) ;;
 	*) help_sections=( $help_switch ) ;;	
esac

for help_section in "${help_sections[@]}"
	do
		source "$this_script"_print_help_"$help_section".sh
done

cat << EOF
======================================================================
EOF
