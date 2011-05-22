
FC=$(FC90)

fname = sdrive.o lbfgs.o 

main : $(fname)
	$(FC) $(fname) $(extras) -o sdrive

.f.o :; $(FC) -c $*.f

