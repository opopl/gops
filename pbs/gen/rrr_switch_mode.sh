
bare_inner_title="T=$temp r=$radius f=$force nrg=$nrg n=$nsteps date= `date_dm_hm`"
od=out-T_$temp-N_$nsteps-nrg_$nrg-r_$radius-f_$force
ofs=( r.log gmin_out ac_ratios markov.dat energy.dat ef.dat gm.xyz )

case "$rrr_mode" in
	"rr") 
		echo "Command: rr 0 -nsteps $nsteps -q $queue -temp $temp -force $force -od $od" >> rrr.log	
		rr 0 -nsteps $nsteps -q $queue -temp $temp -force $force -od $od	
		sleep 3
	;;
	"cp")
		for of in "${ofs[@]}"
			do
				#echo "@@ from dir $od" >> out/$of
				case "$of" in
				  	"markov.dat" | "energy.dat" )
					cat $od/$of >> out/${of%.dat}.$od.dat ;;
					*) 
					cat $od/$of >> out/${of%.dat}.$od 
					cat $od/$of >> out/$of 
					;;
				esac
		done
	;;
	"pem") # {{{

	output_dir=out
	plot_dir=plots
	df_insert=".$od"
	view_plot=0
	inner_title="$bare_inner_title"

	source plot_energy_markov_dat.sh
                                                                                 
	cd $plot_dir
	pdf_add enq.pdf enq$df_insert.pdf
	pdf_add maq.pdf maq$df_insert.pdf
	cd -
	;;
	# }}}
	"pemf")  # {{{
				cd $plot_dir
				pdf_merge em.pdf enq.pdf maq.pdf  
				scp em.pdf $leonov:~/plots/em-$main_mode-`date_dm_hm`.pdf
				scp em.pdf $leonov:~/plots/em-last.pdf
				cd -
	;;
	# }}}
esac
