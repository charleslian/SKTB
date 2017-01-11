///\brief Derive Fitting_band class from Fitting class. 
///\details class to define a specific method to fit parameters using energy band data. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 


#if !defined(_FITTING_BAND_H)
#define _FITTING_BAND_H

#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 
#include "Fitting.h"
#include "Hamiltonian_regular.h"
#include "Parameters.h"
#include "Eigenvalue.h"

class Fitting_band : public Fitting
{
public:
	///\brief constructor for band fitting. 
	///\param _fitting_hamiltonian the hamiltonian to calculate energy band. 
	///\param _parameters_band the parameters to calculate energy band. 
	///\param _eigenvalue_target the target energy band. 
	///\param _lowest_band the lowest band for fitting. 
	///\param _highest_band the highest band for fitting. 
	///\param _criteria the criteria of fitting end condition. 
	Fitting_band(Hamiltonian_regular &_fitting_hamiltonian, Parameters &_parameters_band, Eigenvalue &_eigenvalue_target, int &_lowest_band, int &_highest_band, double &_criteria);
	
	~Fitting_band();///<destructor. 
private:
	int lowest_band;///<lowest band for fitting. 
	int highest_band;///<highest band for fitting. 
	Hamiltonian_regular fitting_hamiltonian;///<the hamiltonian to calculate energy band.  

	///\brief set the vertex despite of control. 
	///\param _index_of_vertex is the index of vertex.  
	///\return the vertex without control. 
	Eigen::MatrixXd set_simplex_vertex_without_control(int _index_of_vertex);

	///\brief function to be optimized. 
	///\return the function value at current parameters. 
	double fitting_function();
};

#endif  //_FITTING_BAND_H
