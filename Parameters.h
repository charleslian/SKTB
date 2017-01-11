///\brief define Parameters. 
///\details class to supply some interfaces for the calculations followed. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 

#if !defined(_PARAMETERS_H)
#define _PARAMETERS_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include<algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 

class Parameters
{
public:
	std::string *parameter_name;///<parameter name. 
	Eigen::MatrixXd values;///<parameter values.
	Eigen::MatrixXd fitting_control;///<fitting control(1 fitting, 0 not fitting). 

	int output(char *_filename);///<output values to _filename. 
	
	Parameters();///<constructor. 
	
	~Parameters();///<destructor. 
};

#endif  //_PARAMETERS_H
