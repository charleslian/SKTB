///\brief Derive Atoms_vasp using Atoms class. 
///\details class supply a implementation to read the data in Atoms from VASP. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 


#if !defined(_ATOMS_VASP_H)
#define _ATOMS_VASP_H

#include "Atoms.h"

class Atoms_vasp : public Atoms
{
public:

	///\brief constructor. 
	///\param _nth_neighbor is the highest order of neighbors to calculate. 
	Atoms_vasp();
	
	~Atoms_vasp();///<\brief destructor. 

private:
	///\brief read *atom and unit_cell_basis from POSCAR
	///\return 0 if succeed. 
	int read();

};

#endif  //_ATOMS_VASP_H
