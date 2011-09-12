
cat << EOF	
SCRIPTS
	hr 			display help info about runs
	rr0			make runs on multiple procs for diff. forces 
				output gathered together

	rr OPTION		shortcut for r -r OPTION

	rrr, sub.sh		PBS job submission

	c			clean
	vf			view files
	vo			view output
	
	qs			view combined output from qstat -q, showq
	qa			view my jobs
	qk			kill my jobs

	crg			compute radii of gyration, and write
					them into "rgyr" output directory
	
	cpt			copy different things

	------------
	PLOTTING
	------------	
	vpl			view a plot with ps/pdf viewer
	pl			plot smth
	gp			plot ef, rgyr

	ao			make an archive copy of all output files
	rmo			delete output dir copy with index i
					(see ao.log)
EOF
