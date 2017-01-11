///\brief Derive Hamiltonian_tight_binding class from Hamiltonian class. 
///\details class to implement Hamiltonian class by dividing the hamiltonian by atom block. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 


#if !defined(_HAMILTONIAN_TIGHT_BINDING_H)
#define _HAMILTONIAN_TIGHT_BINDING_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 
#include "Hamiltonian.h"
#include "Atoms.h"
#include "K_points.h"
#include "Parameters.h"
#include "Hamiltonian.h"
#include "Eigenvalue.h"
#include "Control.h"

class Hamiltonian_tight_binding : public Hamiltonian
{
public:
	int read_control(Control & _control);

	///\brief read data from Atoms class. 
	///\param _atoms the Atoms class contains the structure data. 
	///\return 0 if succeed.
	int read_atoms(Atoms &_atoms);

	///\brief read data from K_points class. 
	///\param _k_points the K_points class contains the k-points data. 
	///\return 0 if succeed.
	int read_k_points(K_points &_k_points);

	///\brief read data from Parameters class. 
	///\param _parameters the Parameters class contains the parameters data. 
	///\return 0 if succeed. 
	int read_parameters(Parameters &_parameters);

	///\brief calculate the hamiltonian. 
	///\return 0 if succeed.
	int calculate();

	int calculate_eigenvalue(Eigen::MatrixXcd hamiltonian,int _k_points_index);
	
	int calculate_DOS();

	int init(Atoms &_atoms, K_points &_k_points, Parameters &_parameters, Control &_control);///<constructor. 

	Hamiltonian_tight_binding();///<destructor. 
	~Hamiltonian_tight_binding();///<destructor. 

	int output(char * filename);
	Eigenvalue eigenvalue;

protected:
	Atoms atoms;///<store the data read from _atoms. 
	K_points k_points;///<store the data read from _k_points. 
	Parameters parameters;///<store the data read from _parameters.
	Control control;
	double max;///<the max order of neighbor. 
	int block_size;///<the rows of the block of hamiltonian. 
	int block_num;///<the row index of the block of hamiltonian. 
	///\brief calculate hamiltonian of k. 
	///\param _k_points_index the index of k points. 
	///\return the hamiltonian of k. 
	Eigen::MatrixXcd calculate_hamiltonian_of_k_index(int _k_points_index);

	///\brief split the hamiltonian by blocks; Each block a pair of atom-atom interaction.  
	///\return 0 if succeed. 
	virtual int split_hamiltonian();

	///\brief calculate self-energy term, i.e. the i-i block. 
	///\param _atom atoms::atom[i]. 
	///\param _k_points_index is the k point. 
	///\return the i-i block hamiltonian of k. 
	virtual Eigen::MatrixXcd calculate_self_energy(Atoms::Atom _atom, int _k_points_index);

	///\brief calculate calculate hooping term, i.e. i-j block. 
	///\param _atom1 atoms::atom[i]. 
	///\param _atom2 atoms::atom[j]. 
	///\param _k_points_index is the k point. 
	///\return the i-j block hamiltonian of k. 
	virtual Eigen::MatrixXcd calculate_hooping(Atoms::Atom _atom1, Atoms::Atom _atom2, int _k_points_index);

	double fermi_smear(double energy,double ISMEAR);
};

#endif  //_HAMILTONIAN_TIGHT_BINDING_H
