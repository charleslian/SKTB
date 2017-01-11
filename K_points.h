///\brief Define K_points class. 
///\details class to supply some interfaces for the calculations followed. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 

#if !defined(_K_POINTS_H)
#define _K_POINTS_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 

class K_points
{
public:
	Eigen::RowVectorXd band_x_coordinates;///<the the x coordinates(distance) of the energy band. 
	Eigen::MatrixXd k_points_coordinates;///<the the k points coordinates. 
	int k_points_num;///<the number of k-points. 
	
	int output();///<output the band_x_coordinates and k_points_coordinates to k_coordinates.dat. 
	
	K_points();///<\brief constructor. 
	
	~K_points();///<\brief destructor. 
};
#endif  //_K_POINTS_H
