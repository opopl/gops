
if [ $reset_pf -eq 1 ]; then
	plot_file=$plot_dir/$arg.$ext
else
  	plot_file=$plot_dir/$pf.$ext
fi

gp_terminal="postscript eps enhanced"

case "$plot_type" in # {{{
  	"plot") # {{{ 
data_file=$output_dir/$ifile.dat 	
gm_en=-0.2936531215
e0=$gm_en

echo "$plot_file"

rm -rf $plot_file

case "$arg" in
  	"ef") # energy vs force {{{
		case "$gpi" in
	  	"0") # {{{
cat > $arg.gp << EOF

set pointsize 2
unset key
set term $gp_terminal 
set output "$plot_file" 

set xlabel 'f'
set ylabel 'E(f)'

plot "$data_file"  using 1:6 with points pt 4 ps 2 

EOF
	;;
	# }}}
	"logxy") # {{{ 
		fh=` head -1 "$data_file" | sed '/#/d' `

		#temp=echo $fh | awk '{ print } '
		
		fmin=` awk 'NR==is { print $colx }' colx=$colx is=$is $data_file `
		emin=` awk 'NR==is { print $coly }' coly=$coly is=$is $data_file `

		fmax=` awk 'END { print $colx }'  colx=$colx  $data_file `
		emax=` awk 'END { print $coly }'  coly=$coly $data_file `

		#emin=$gm_en

		#df=$( echo "($fmin)" | bc )
		#df=$( echo "1.0" | bc )
		de=$( echo "($gm_en)" | bc )

		echo "Start energy: emin=$emin" 
		echo "Start force:  fmin=$fmin"
		
		echo "End energy: emax=$emax" 
		echo "End force:  fmax=$fmax"

		echo "df=$df"
		echo "de=$den"
		
		cat $data_file | sed '/#/d' | awk ' 
	       		{ f=$colx; printf "%s ",f
			  for (i=coly;i<=NF;i++) 
			  {ef=($i-e0)/de; printf "%s ", ef;} 
			  printf "\n"; } ' colx=$colx coly=$coly e0=$gm_en de=$de > n

		cat n | sed '/#/d' | awk ' 
		 	NR==is,0 { for (i=1;i<=2;i++) {ef=$i; printf "%s ", log(ef)/log(10);} 
			printf "\n"; } ' is=$is > n1 

		data_file=n1

cat >& fit.gp << EOF
f(t)=a*t+b
fit f(x) "$data_file" using 1:2 via a,b 
EOF

gnuplot fit.gp >& fit.tmp 

awk '/Final set of parameters/,/correlation matrix/ { print $0 }' fit.tmp >& fit.tmp.1
mv fit.tmp.1 fit.tmp

a_fit=` awk '/^a[ ]*=/ { print $3 }'  fit.tmp `
a_err=` awk '/^a[ ]*=/ { print $5 }'  fit.tmp  `

b_fit=` awk '/^b[ ]*=/ { print $3 }'  fit.tmp  `
b_err=` awk '/^b[ ]*=/ { print $5 }'  fit.tmp  `

a_fit=` printf "%3.2f\n" $a_fit `
b_fit=` printf "%3.2f\n" $b_fit `

a_err=` printf "%3.2f\n" $a_err `
b_err=` printf "%3.2f\n" $b_err ` 

cat > $arg.gp << EOF

set pointsize 2
set key left bottom

set term $gp_terminal "$gp_font"
set output "$plot_file" 

set title "Global minimum energy vs applied pulling force"

set ylabel 'log_{10}{/Symbol e}'          
set xlabel 'log_{10}f'
set ylabel 'log_{10}(E-E_0)/E_0'

EOF

labels=( "Data File: $ifile" "f_{min}=$fmin" "f_{max}=$fmax" "E_0=$gm_en" "E_{min}=$emin" "E_{max}=$emax" "Curve fit: \n a=$a_fit, da=$a_err \n b=$b_fit, db=$b_err" )

let i=0
for label in "${labels[@]}"
	do
	  dy=$( echo "-$i*$lshift" | bc )
	  position="$lk $lx, $lk $ly+$dy"
	  echo "set label \"$label\" at $position" >> $arg.gp
	  i=$(($i+1))
done

cat >> $arg.gp << EOF

f(x)=$a_fit*x+$b_fit
set size $gp_xsize,$gp_ysize
plot "$data_file" using 1:2  title "Data points" with points pt 4 ps 1, f(x) title "Linear least-squares fit, f(x)=ax+b"

EOF
	# }}}
	;;
esac

gnuplot $arg.gp
cp $plot_file $rpdir

;;
# }}}
	"s_f0")  # energy.dat, markov. dat {{{ 
	df_insert=""
	source plot_energy_markov_dat.sh
;;
# }}}
esac ;;
# }}}
	"xyz") # {{{
		data_file=$output_dir/$arg.xyz
		pymol -c $arg.pml $data_file 
	;;
	# }}}
esac # }}}
