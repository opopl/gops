
gp_terminal="postscript enhanced portrait"
gp_terminal="postscript enhanced landscape"

data_files=( energy markov )

for df in ${data_files[@]} 
	do
	  data_file=$output_dir/$df$df_insert.dat
	  plot_file=$plot_dir/$df.ps
	  case "$df" in
		"energy") 
			plot_title="Energy vs Quench Number $add_plot_title" 
			y_label="Energy"
			;;
		"markov") 
			plot_title="Markov Energy vs Quench Number $add_plot_title" 
			y_label="Markov Energy"
			;;
	  esac
	  
it_x=0.05
it_y=0.95

cat > em.gp << EOF

set style line 1 lt 2 lw 3
set key box linestyle 1

set term $gp_terminal 
set output "$plot_file" 
set title "$plot_title"
set xlabel "Quench Number"
set ylabel "$y_label"

set label "$inner_title"  at graph  $it_x, graph  $it_y

plot "$data_file" using 1:2 $o notitle

EOF
gnuplot em.gp

done
cd $plot_dir

#gs -dNOPAUSE -sDEVICE=pswrite -sOutputFile=$arg.ps energy.ps markov.ps -c quit

ps2pdf energy.ps enq$df_insert.pdf
ps2pdf markov.ps maq$df_insert.pdf

#pdf_merge $arg.pdf energy.pdf markov.pdf

cd -
plot_file=$plot_dir/$arg.pdf
gv_opts="-landscape"

