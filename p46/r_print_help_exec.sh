cat << EOF
EXECUTION REGIMES
	-h, -h*				print (this) full help message (on different help topics)
						* can be=
					       		empty string (for everything)	
							r (for execution regimes), 		
	-pef,	--plot_ef		plot energy vs force, only, using input file ef.dat
	-pgm,	--plot_gm_xyz		plot 3D pics of all gm minimuma in gm.xyz
					into a movie
	
	-r regime_name			specific execution regime. For more info, 
						supply -hr help option
EOF

