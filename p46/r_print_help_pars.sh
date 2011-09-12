
cat << EOF
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
EOF
