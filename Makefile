target=Atoms.o Atoms_vasp.o Control.o Eigenvalue.o Eigenvalue_ab_initio.o Fitting.o Fitting_band.o Hamiltonian.o Hamiltonian_tight_binding.o Hamiltonian_regular.o  Hamiltonian_SOC.o K_points.o K_points_grid.o K_points_line.o Parameters.o Parameters_band.o Slater_Koster.o

fopenmp=-fopenmp

slater-koster:$(target)
	g++ -o slater-koster $(fopenmp) $(target)

Atoms.o:Atoms.cpp Atoms.h
	g++ -c -O3 $(fopenmp) Atoms.cpp

Atoms_vasp.o:Atoms_vasp.cpp Atoms_vasp.h
	g++ -c -O3 $(fopenmp) Atoms_vasp.cpp

Control.o:Control.cpp Control.h
	g++ -c -O3 $(fopenmp) Control.cpp

Eigenvalue.o:Eigenvalue.cpp Eigenvalue.h
	g++ -c -O3 $(fopenmp) Eigenvalue.cpp

Eigenvalue_ab_initio.o:Eigenvalue_ab_initio.cpp Eigenvalue_ab_initio.h
	g++ -c -O3 $(fopenmp) Eigenvalue_ab_initio.cpp


Fitting.o:Fitting.cpp Fitting.h
	g++ -c -O3 $(fopenmp) Fitting.cpp

Fitting_band.o:Fitting_band.cpp Fitting_band.h
	g++ -c -O3 $(fopenmp) Fitting_band.cpp

Hamiltonian.o:Hamiltonian.cpp Hamiltonian.h
	g++ -c -O3 $(fopenmp) Hamiltonian.cpp

Hamiltonian_tight_binding.o:Hamiltonian_tight_binding.cpp Hamiltonian_tight_binding.h
	g++ -c -O3 $(fopenmp) Hamiltonian_tight_binding.cpp

Hamiltonian_regular.o:Hamiltonian_regular.cpp Hamiltonian_regular.h
	g++ -c -O3 $(fopenmp) Hamiltonian_regular.cpp

Hamiltonian_SOC.o:Hamiltonian_SOC.cpp Hamiltonian_SOC.h
	g++ -c -O3 $(fopenmp) Hamiltonian_SOC.cpp

K_points.o:K_points.cpp K_points.h
	g++ -c -O3 $(fopenmp) K_points.cpp

K_points_line.o:K_points_line.cpp K_points_line.h
	g++ -c -O3 $(fopenmp) K_points_line.cpp

K_points_grid.o:K_points_grid.cpp K_points_grid.h
	g++ -c -O3 $(fopenmp) K_points_grid.cpp

Parameters.o:Parameters.cpp Parameters.h
	g++ -c -O3 $(fopenmp) Parameters.cpp

Parameters_band.o:Parameters_band.cpp Parameters_band.h
	g++ -c -O3 $(fopenmp) Parameters_band.cpp

Slater_Koster.o:Slater_Koster.cpp
	g++ -c -O3 $(fopenmp) Slater_Koster.cpp 

clean: 
	rm  $(target)

clean-all:
	rm $(target) Slater-Koster
