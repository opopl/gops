
cat << EOF
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
	-sys system			specify model system.
						Values for variable system:
							P46 (wild-type BLN model)
							G46 (Go-like BLN model)
						Default value is G46.
EOF
