# Compiler choice
FC = gfortran
MPIFC = gfortran
IFC = ifort

MLIBS =

#################

EXES = tpcf.exe ccf.exe

all: $(EXES)

tpcf.exe: tpcf.f90
	$(FC) -fopenmp tpcf.f90 -o tpcf.exe $(MLIBS)

ccf.exe: ccf.f90
	$(FC) -fopenmp ccf.f90 -o ccf.exe $(MLIBS)


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
