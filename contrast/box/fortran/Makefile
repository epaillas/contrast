# Compiler choice
FC = gfortran
MPIFC = gfortran

MLIBS =

#################

EXES = tpcf.exe tpcf_2d.exe tpcf_ap.exe

all: $(EXES)

tpcf.exe: tpcf.f90
	$(FC) -fopenmp tpcf.f90 -o tpcf.exe $(MLIBS)

tpcf_2d.exe: tpcf_2d.o
	$(FC) tpcf_2d.o -o tpcf_2d.exe $(MLIBS)

tpcf_ap.exe: tpcf_ap.o
	$(FC) tpcf_ap.o -o tpcf_ap.exe $(MLIBS)
	

tpcf_2d.o: tpcf_2d.f90
	$(FC) -c tpcf_2d.f90

tpcf_ap.o: tpcf_ap.f90
	$(FC) -c tpcf_ap.f90


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
