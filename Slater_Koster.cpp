// Slater_Koster.cpp
#define PI 3.141592653
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <omp.h>
#include <time.h>

#include "Eigen/Core"
#include "Control.h"
#include "Atoms_vasp.h"
#include "K_points_line.h"
#include "K_points_grid.h"
#include "Parameters_band.h"
#include "Hamiltonian_regular.h"
#include "Hamiltonian_SOC.h"
#include "Eigenvalue_ab_initio.h"
#include "Fitting_band.h"

int main(int argc, char* argv[])
{
	std::cout<<"The number of threads is "<<Eigen::nbThreads()<<std::endl;
	std::cout<<"The number of processes is "<<omp_get_num_procs()<<std::endl;
	clock_t start=clock();

	Eigen::initParallel();
	Control control("control");//the control tags. 
	Atoms_vasp atoms;//atoms read from VASP. 
	K_points_line linear_k_points;
	linear_k_points.init();//k-points read from VASP.
	
	Eigenvalue_ab_initio ab_initio_eigenvalue;
	ab_initio_eigenvalue.calculate_eigenvalue();
	ab_initio_eigenvalue.output("ab_initio_energy_band.dat");

	Parameters_band parameters(atoms);
	if (control.fitting_needed)
	{
		Parameters_band parameters(atoms);
		K_points_line linear_k_points_fit;
		linear_k_points.select_k_points(linear_k_points_fit,control.num_of_fitting_point);
		Eigenvalue ab_initio_eigenvalue_selected,target;
		ab_initio_eigenvalue.select_k_points(ab_initio_eigenvalue_selected,control.num_of_fitting_point);
		ab_initio_eigenvalue_selected.select_bands(target,control.fitting_band[0],control.fitting_band[1]);
		Hamiltonian_regular hamiltonian(atoms,linear_k_points_fit,parameters,control);
		Fitting_band fitting(hamiltonian,parameters,target,control.fitting_band[0],control.fitting_band[1], control.criteria);
	}
	if (control.band_calculation_needed)
	{
		Hamiltonian_regular hamiltonian(atoms,linear_k_points,parameters,control);//calculate the hamiltonian. 
		hamiltonian.output("tight_binding_energy_band.dat");
		hamiltonian.eigenvalue.calculate_band_gap_and_Vf(atoms.total_electron_num,control.smear_sigma);
	}
	if (control.SOC_calculation_needed)
	{
		Hamiltonian_SOC hamiltonian_SOC(atoms,linear_k_points,parameters,control);//calculate the hamiltonian. 
		hamiltonian_SOC.output("tight_binding_energy_band_SOC.dat");
	}
	if (control.DOS_calculation_needed)
	{
		K_points_grid k_points_grid(control.k_mesh_num);
		Hamiltonian_regular hamiltonian(atoms,k_points_grid,parameters,control);
		hamiltonian.calculate_DOS();
	}
	clock_t finish=clock();
	double totaltime=(double)(finish-start)/(double)CLOCKS_PER_SEC;
	std::cout<<"The running time is "<<totaltime<<"seconds"<<std::endl;
} 
