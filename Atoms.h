///\brief Define Atoms class. 
///\details class supply some interfaces for the calculations followed. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 

#if !defined(_ATOMS_H)
#define _ATOMS_H

#include <list>
#include "Eigen/Dense"

class Atoms
{
public:
	///\brief Define the struct Element. 
	struct Element
	{
		int index;///<\brief The index of the element. 
		std::string name;///<\brief The name of the element. 
		int electron_number;
	};

	///\brief Define the struct Atom. 
	struct Atom
	{
		Element element;///<\brief The element of this atom. 
		int index;///<\brief The index of this atom in the unit cell. 
		Eigen::RowVector3d coordinates;///<\brief The coordinates (x,y,z) of this atom. 
		//std::list<Atom> neighbors;///<\brief The neighbor list, neighbors[i] store the ith order neighbor list of this atom. 
	};

	struct Neighbor
	{
		int index_i;
		int index_j;
		Eigen::RowVector3d coordinates1;
		Eigen::RowVector3d coordinates2;
	};
	std::list<Neighbor> neighbors;
	Atom *atom;///<\brief The array of atom in the unit cell. 
	Element *element;///<\brief The array of element in the unit cell. 
	int total_atom_num;///<\brief The total number of the atom. 
	int element_species;///<\brief the number of element species. 
	int total_electron_num;
	Eigen::Matrix3d unit_cell_basis;///<\brief unit cell basis. 
	Eigen::Matrix3d reciprocal_basis;///<\brief reciprocal basis. 
	Eigen::MatrixXd cutoff_distance_between_elements;///<\brief neighbor_distance(i,j) store the kth order neighbor list between element i and elment j. 

	///\brief calculate the neighbors of all the atoms in the unit cell. 
	///\attention must be called after unit_cell_basis, element_species, total_atom_num and *atom are initialized. 
    ///\param _nth_neighbor is the highest order of neighbors to calculate. 
    ///\return 0 if succeed. 
	int calculate_neighbor();

	///\brief calculate reciprocal basis and write to REPVEC
	///\return 0 if succeed. 
	int calculate_reciprocal_basis();

	Atoms();///<\brief constructor. 

	~Atoms();///<\brief destructor. 
};

#endif  ///_ATOMS_H
