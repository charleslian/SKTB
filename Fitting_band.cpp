#include "Fitting_band.h"
#include "Hamiltonian_regular.h"
#include "Parameters.h"
#include "Eigenvalue.h"
#include <time.h>

Fitting_band::~Fitting_band() 
{

}

Fitting_band::Fitting_band(Hamiltonian_regular &_fitting_hamiltonian, Parameters &_parameters_band, Eigenvalue &_eigenvalue_target, int &_lowest_band, int &_highest_band, double &_criteria)
{
	//generate seed using system time. 
	srand(unsigned(time(0)));

	//read data to calculate fitting function. 
	fitting_hamiltonian=_fitting_hamiltonian;
	parameters=_parameters_band;
	lowest_band=_lowest_band;
	highest_band=_highest_band;
	target=_eigenvalue_target.eigenvalues;

	//fitting. 
	fitting(_criteria);
	
	//output the fitting results. 
	_parameters_band=parameters;
	parameters.output("NEWPARAM");
}

Eigen::MatrixXd Fitting_band::set_simplex_vertex_without_control(int _index_of_vertex) 
{
	//build the simplex_vertex_without_control.
	int NX=parameters.values.rows();
	int NY=parameters.values.cols();
	int N=NX*NY;
	Eigen::MatrixXd simplex_vertex_without_control(NX,NY);

	//the first vertex is built from original parameters. 
	if(_index_of_vertex==0)
	{
		simplex_vertex_without_control=parameters.values;
		return simplex_vertex_without_control;
	}

	//the vertex is created by adding a random array to the first vertex.  
	Eigen::MatrixXd delta(NX,NY);
	delta.setRandom();
	delta.col(0)*=10;
	simplex_vertex_without_control=parameters.values+delta;
	return simplex_vertex_without_control;
}

double Fitting_band::fitting_function() 
{
	//recalculate the hamiltonian using current parameters.
	fitting_hamiltonian.read_parameters(parameters);
	fitting_hamiltonian.calculate();
	
	//calculate the eigenvalues of the fitting_hamiltonian. 
	Eigenvalue eigenvalue_tight_binding_selected;

	//select certain bands for fitting. 
	fitting_hamiltonian.eigenvalue.select_bands(eigenvalue_tight_binding_selected,lowest_band,highest_band);
	current=eigenvalue_tight_binding_selected.eigenvalues;

	//the fitting function value is the average distance between current and target. 
	Eigen::ArrayXXd fitting_matrix(current.rows(),current.cols());
	fitting_matrix=(current-target).array().abs();
	double fitting_function_value=fitting_matrix.sum()/fitting_matrix.size();
	return fitting_function_value;
}

