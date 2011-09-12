
#gp_opts="with linespoints title \"$title\""
gp_opts="with linespoints title \"$title\""
gp_opts="with linespoints "
gp_opts="with points "

case "$prefix" in
  	"ef") xc=1 ; yc=6 ; xl="f" ; yl="E"
	;;
  	"rgyr") xc=2 ; yc=3 ; xl="f" ; yl="R_g"
       	;;
esac


cat >& $a.gp << EOF

set term postscript enhanced
set output "$output_file"
set title "$title"

set xlabel "$xl"
set ylabel "$yl"

plot "$data_file" using $xc:$yc $gp_opts

EOF

gnuplot $a.gp
