
eo "will try to plot energies vs force..." 
pref="plot.energies"

data=$output_dir/ef.tex

eo "Input data file: $data" 

o="with lines"
o="with linespoints"

mkdir -p plots
plot_file=plots/ef.ps

rm -rf $plot_file

cat > ef.gp << EOF

set term postscript enhanced  
set output "$plot_file" 
set xlabel 'Force'          
set ylabel 'Energy'  
set label "N=$nsteps" at 2,0
#set multiplot  
plot "$data" using 1:2 $o, "$data" using 1:3 $o, "$data" using 1:4 $o 
#set nomultiplot 

EOF

#cat ef.gp | gnuplot >> $logf
gnuplot ef.gp

if [ -f "$plot_file" ]; then
		f_size=$(stat -c%s "$plot_file")
  		eo "Output plot file: $plot_file ($f_size)"
      else
		eo "Didn't produce $plot_file"
fi	

#ps2pdf ef.ps ef.pdf

