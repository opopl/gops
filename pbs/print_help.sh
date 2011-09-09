

case "$help_switch" in
  	"all")
cat << EOF

======================================================================
FILENAME: r
PURPOSE: bash executable script which submits jobs to the PBS queues
USAGE: r [ OPTIONS ]

==========================
OPTIONS
==========================
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

PBS 

	-q, --pbs_queue_name string	queue name
	
	--pbs_job_name 	 string		job name; default is the directory name 
					in which script is located

	-wt, --pbs_wall_time	string	specify PBS walltime
						  "day"	24:00:00

	============================
	include inside ...> in r.log:
	============================

	--pref_pbs_job_id		pbs_job_id
	--pref_pbs_job_name		pbs_job_name, specified by the --pbs_job_name switch
	--pref_job_number		job_number, specified by the -j switch
GMIN
	--print_gmin_warnings		print warnings from GMIN output file GMIN_out

	--remove_dn
	--not_remove_dn

	-bs, --basin-sampling		specify a basin-sampling run. Keywords in the 'data'
						file:
							HISTOGRAM
							TEMPERATURE 0
							BINSTRUCTURES
							TETHER

JOBS
	-j job_number		assign to this job a specific job_number (integer => 1)
					e.g., r -j 2 means all output files will be appended
					job_number=2 at the end, e.g., ef.dat => ef.dat.2 etc.
OUTPUT OPTIONS
	varying output filenames:

	--append_job_number	append job number, only makes sense if the option -j is 
       					specified	
					set internal parameter append
	--n-append_job_number	don't append job number
	--append_job_name	append job name, specified by the switch --pbs_job_name 

	-od, --output_dir_name	specify the output directory name (default is out)
				
TASKS

	-t-gmgof		run GMIN for Go model with different pulling forces 
					Required keywords in GMIN data file:
						PULL 
						G46
	-t-gmwtf		as -t-gmgof, except for WT model
						in data, keyword changed G46 -> P46
PARAMETERS

	-radius			container radius
	-nrg			number of starting random geometries
	-seed			random number seed
	-nsteps int		number of MC steps

	FOR PULLING:

	--adj_target i		adjust the target energy for each force
					This makes sense as long as do_force_loop=1
					Possible values of i are:
						0 	don't adjust
						1 	target_energy=gm_energy

        -if [ f1 f2 f3 ]	the array of forces is obtained by incrementing
					from f1 to f2, with increment f3
					set internal parameter get_force="increment"
	--no_force_loop		don't loop over forces, for the force use the value
					specified by the option below,
					set internal parameter do_force_loop=0
	--force_loop		loop over forces,
					set internal parameter do_force_loop=1

	-force float		specify a pulling force,
       					set internal parameter get_force="single"

	-af 			forces are specified by an array, tbc
	--modify-te			modify the target energy, specified by keyword TARGET.
       					Values are:
						0 - don't modify
						1 - divide by 1+F, with F applied force

	OTHER

	-te	float		specify target energy. In data file, it corresponds to the keyword
					TARGET	
					default value: -0.2936531215

	-temp float			temperature, default value 0.03

	-n_en int			number of lowest energies to extract 
					from file lowest, default 4
	-pstart int		starting pulling point, default 1
	-pend 	int		end pulling point,	default 46

	-ap,	--append_files	
				append output files, don't overwrite
					sets append_files=true
					(default: append_files=false)
	-nap,	--not_append_files	
				don't append output files
					append_files=false

OTHER SCRIPTS:
	
	hr 			display help info about runs
	rr0			make runs on multiple procs for diff. forces 
				output gathered together

	rr OPTION		shortcut for r -r OPTION

	c			clean
	vf			view files
	vo			view output
	
	qs			view combined output from qstat -q, showq
	qa			view my jobs
	qk			kill my jobs
	
	vpl			view a plot with ps/pdf viewer
	pl			plot smth

	ao			make an archive copy of all output files
	rmo			delete output dir copy with index i
					(see ao.log)

AUTHOR: O. Poplavskyy
======================================================================

EOF
;;
"runs")
	echo "======================================================================"
	echo "	RUNS:"
	let i=0
	for rn in "${runs[@]}"
		do
		  echo "$i	$rn"
		  i=$(($i+1))
	done
	echo "======================================================================"
;;
esac

