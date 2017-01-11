///\brief Derive K_points_grid class from K_points class. 
///\details class to implement the K_points class using Monkhorst-Pack method. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 

#if !defined(_K_POINTS_GRID_H)
#define _K_POINTS_GRID_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include<algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 
#include "K_points.h"

class K_points_grid : public K_points
{
public:
	///\brief constructor. 
	///\param _k_mesh_num the mesh-grid number in three direction. 
	K_points_grid(Eigen::RowVector3d &_k_mesh_num);

	K_points_grid();///<\brief constructor. 

	~K_points_grid();///<\brief destructor. 

private:

	Eigen::Matrix3d reciprocal_basis;///<\brief reciprocal basis. 

	int read_reciprocal_basis();///<\brief read reciprocal basis from REPVEC. 

	///\brief calculate k grid. 
	///\param _k_mesh_num the mesh-grid number in three direction. 
	int calculate_k_grid(Eigen::RowVector3d &_k_mesh_num);
};

#endif  //_K_POINTS_GRID_H
