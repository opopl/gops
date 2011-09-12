F90=$(FC90)
F90FLAGS = -O

all: disconnectionDPS manipulate

disconnectionDPS: disconnectionDPS.f90
	$(F90) $(F90FLAGS) -o $@ $<
#	cp $@ ddps

manipulate: manipulate.f90
	$(F90) $(F90FLAGS) -o $@ $<
#	cp $@ mnp

clean:
	rm -f disconnectionDPS manipulate *.o *.mod
