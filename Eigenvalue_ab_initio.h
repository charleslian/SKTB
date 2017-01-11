///\brief Derive Eigenvalue_ab_initio class from Eigenvalue class. 
///\details class to implement the eigenvalue class by reading data from ab-initio results. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 

#if !defined(_EIGENVALUE_AB_INITIO_H)
#define _EIGENVALUE_AB_INITIO_H

#include "Eigenvalue.h"

class Eigenvalue_ab_initio : public Eigenvalue
{
public:
	///\brief read the eigenvalue from EIGENVAL and DOSCAR. 
	///\return 0 if succeed. 
	int calculate_eigenvalue();
	
	Eigenvalue_ab_initio();///<\brief constructor. 

	~Eigenvalue_ab_initio();///<\brief destructor. 
};

#endif  //_EIGENVALUE_AB_INITIO_H
