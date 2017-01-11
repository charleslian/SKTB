///\brief Derive Parameters_band class from Parameters class. 
///\details class to implement the Parameters class by reading PARAM file for band calculation. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 

#if !defined(_PARAMETERS_BAND_H)
#define _PARAMETERS_BAND_H

#include "Parameters.h"
#include "Atoms.h"

class Parameters_band : public Parameters
{
public:
	///\brief constructor. 
	///\param _atoms supply the element species and element names.
	Parameters_band(Atoms &_atoms);

	Parameters_band();///<constructor. 

	~Parameters_band();///<destructor. 
};

#endif  //_PARAMETERS_BAND_H
