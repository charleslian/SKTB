///\brief define Fitting class. 
///\details class to fit parameters using the Nelder-Mead method. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 


#if !defined(_FITTING_H)
#define _FITTING_H

#include "Eigen/Dense"
#include "Parameters.h"

class Fitting
{
public:
	Eigen::MatrixXd current;///<current function value. 
	Eigen::MatrixXd target;///<target function value. 
	Parameters parameters;///<current parameters. 

	///\brief fitting parameters using the Nelder-Mead method. 
	///\param _criteria the criteria of ending the fitting.  
	///\return 0 if succeed. 
	int fitting(double _criteria);

	Fitting();///<constructor. 

	~Fitting();///<destructor. 
protected:
	Eigen::MatrixXd simplex;///<current function value. 

	///\brief set the vertex despite of control. Virtual function which should be overridden in a derived class.
	///\param _index_of_vertex is the index of vertex.  
	///\return the vertex without control. 
	virtual Eigen::MatrixXd set_simplex_vertex_without_control(int _index_of_vertex);

	///\brief create simplex from the vertex without control, the parameters with control=0 is omitted.   
	///\return 0 if succeed. 
	int create_simplex();

	///\brief convert the vertex of the simplex to the parameters. 
	///\param _try_simplex_vertex is the vertex to be converted.  
	///\return 0 if succeed. 
	int simplex_to_parameters(Eigen::RowVectorXd _try_simplex_vertex);

	///\brief function to be optimized. Virtual function which should be overridden in a derived class.
	///\return the function value at current parameters. 
	virtual double fitting_function();
};
#endif  ///_FITTING_H
