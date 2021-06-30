# Compiler choice
FC = gfortran
MPIFC = gfortran

MLIBS =

#################

EXES = tpcf_monopole.exe tpcf_rmu.exe mean_radial_velocity_monopole.exe std_los_velocity_rmu.exe tpcf_2d.exe

all: $(EXES)

tpcf_monopole.exe: tpcf_monopole.f90
	$(FC) -fopenmp procedures.f90 tpcf_monopole.f90 -o tpcf_monopole.exe $(MLIBS)

tpcf_rmu.exe: tpcf_rmu.f90
	$(FC) -fopenmp procedures.f90 tpcf_rmu.f90 -o tpcf_rmu.exe $(MLIBS)

mean_radial_velocity_monopole.exe: mean_radial_velocity_monopole.f90
	$(FC) -fopenmp procedures.f90 mean_radial_velocity_monopole.f90 -o mean_radial_velocity_monopole.exe $(MLIBS)

std_los_velocity_rmu.exe: std_los_velocity_rmu.f90
	$(FC) -fopenmp procedures.f90 std_los_velocity_rmu.f90 -o std_los_velocity_rmu.exe $(MLIBS)

tpcf_2d.exe: tpcf_2d.f90
	$(FC) -fopenmp procedures.f90 tpcf_2d.f90 -o tpcf_2d.exe $(MLIBS)


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
