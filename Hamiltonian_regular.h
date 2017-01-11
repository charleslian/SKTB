///\brief Derive Hamiltonian_regular class from Hamiltonian_tight_binding class. 
///\details class to implement Hamiltonian_tight_binding class by specify the method to split and calculate hamiltonian blocks.  
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 

#if !defined(_HAMILTONIAN_REGULAR_H)
#define _HAMILTONIAN_REGULAR_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include<algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 
#include "Hamiltonian_tight_binding.h"
#include "Atoms.h"
#include "K_points.h"
#include "Parameters.h"

class Hamiltonian_regular : public Hamiltonian_tight_binding
{
public:
	///\brief constructor. 
	///\param _atoms the Atoms class contains the structure data. 
	///\param _k_points the K_points class contains the k-points data. 
	///\param _parameters the Parameters class contains the parameters data. 
	///\param _max the max order of neighbor. 
	Hamiltonian_regular(Atoms &_atoms, K_points &_k_points, Parameters &_parameters, Control &_control);
	
	Hamiltonian_regular();///<constructor. 
	
	~Hamiltonian_regular();///<destructor. 

protected:
	///\brief split the hamiltonian by blocks; Each block a pair of atom-atom interaction.  
	///\return 0 if succeed. 
	int split_hamiltonian();

	///\brief calculate self-energy term, i.e. the i-i block. 
	///\param _atom atoms::atom[i]. 
	///\param _k_points_index is the k point. 
	///\return the i-i block hamiltonian of k. 
	Eigen::MatrixXcd calculate_self_energy(Atoms::Atom _atom, int _k_points_index);
	
	///\brief calculate calculate hooping term, i.e. i-j block. 
	///\param _atom1 atoms::atom[i]. 
	///\param _atom2 atoms::atom[j]. 
	///\param _k_points_index is the k point. 
	///\return the i-j block hamiltonian of k. 
	Eigen::MatrixXcd calculate_hooping(Atoms::Atom _atom1, Atoms::Atom _atom2, int _k_points_index);
};

#endif  //_HAMILTONIAN_REGULAR_H
