///\brief define Hamiltonian class. 
///\details class to supply a interface for the Eigenvalue_tight_binding class to calculate eigenvalues. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 


#if !defined(_HAMILTONIAN_H)
#define _HAMILTONIAN_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include<algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 

class Hamiltonian
{
public:
	Eigen::MatrixXcd hamiltonian;///<hamiltonian.

	int k_points_num;///<the number of k points. 
	
	Hamiltonian();///<constructor. 
	
	~Hamiltonian();///<destructor. 
};

#endif  //_HAMILTONIAN_H
