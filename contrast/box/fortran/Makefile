# Compiler choice
FC = gfortran
MPIFC = gfortran

MLIBS =

#################

EXES = tpcf.exe mean_radial_velocity_vs_r.exe std_los_velocity_vs_rmu.exe tpcf_2d.exe tpcf_ap.exe

all: $(EXES)

tpcf.exe: tpcf.f90
	$(FC) -fopenmp procedures.f90 tpcf.f90 -o tpcf.exe $(MLIBS)

mean_radial_velocity_vs_r.exe: mean_radial_velocity_vs_r.f90
	$(FC) -fopenmp procedures.f90 mean_radial_velocity_vs_r.f90 -o mean_radial_velocity_vs_r.exe $(MLIBS)

std_los_velocity_vs_rmu.exe: std_los_velocity_vs_rmu.f90
	$(FC) -fopenmp procedures.f90 std_los_velocity_vs_rmu.f90 -o std_los_velocity_vs_rmu.exe $(MLIBS)

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
