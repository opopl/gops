
plot_file=$plot_dir/$arg.$ext
gp_terminal="postscript enhanced"

case "$plot_type" in
  	"plot") 
data_file=$output_dir/$arg.dat 	

rm -rf $plot_file

#o="with lines"
#o="with linespoints"
o="with points"

case "$arg" in # {{{
  	"gm") # global minima configurations
	
	;;
  	"ef") # energy vs force {{{
cat > $arg.gp << EOF

set pointsize 2
set term $gp_terminal 
set output "$plot_file" 
set xlabel 'Force'          
#set xrange [0.001:0.005]
set ylabel 'Energy'  
set label "N=$nsteps" at 0,0.5
set title "$plot_title"
# at 0.001,-0.5
#set label "Output dir: $output_dir" at 0.001,-0.55
plot "$data_file" using 1:2 $o, "$data_file" using 1:3 $o, "$data_file" using 1:4 $o 

EOF
;;
# }}}
	"s_f0")  # energy.dat, markov. dat {{{ 
	df_insert=""
	source plot_energy_markov_dat.sh
;;
# }}}
esac 
# }}}
;;
	"xyz")
		data_file=$output_dir/$arg.xyz
		pymol -c $arg.pml $data_file 
	;;
esac

