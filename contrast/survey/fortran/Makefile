# Compiler choice
FC = gfortran
MPIFC = gfortran
IFC = ifort

MLIBS =

#################

EXES = tpcf.exe tpccf.exe

all: $(EXES)

tpcf.exe: tpcf.f90
	$(FC) -fopenmp procedures.f90 tpcf.f90 -o tpcf.exe $(MLIBS)

tpccf.exe: tpccf.f90
	$(FC) -fopenmp procedures.f90 tpccf.f90 -o tpccf.exe $(MLIBS)


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
