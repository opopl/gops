# > ifort -V
# Intel(R) Fortran Compiler for applications running on Intel(R) 64, Version 10.1    Build 20081024 Package ID: l_fc_p_10.1.021
# Copyright (C) 1985-2008 Intel Corporation.  All rights reserved.

vpath %.f90 src

OBJ := bowman3b.o hb3b_coef.o hb3b.o p12body.o pes_shell.o bowmanwater.o 

# Targets:
.PHONY : clean 

../libbowman.a: $(OBJ)
	ar sr $@ $^ ./libs/libpes.a ./libs/libpesd.a

%.o : %.f90
	ifort -r8 -c -O -u -I./mod $< 

clean :
	rm -f *.o *.mod 
