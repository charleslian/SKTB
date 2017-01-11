///\brief Derive K_points_line class from K_points class. 
///\details class to implement the K_points class by reading Ab-initio data. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 


#if !defined(_K_POINTS_LINE_H)
#define _K_POINTS_LINE_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include "Eigen/Dense"
#include "Eigen/Geometry" 
#include "Eigen/Eigenvalues" 
#include "K_points.h"

class K_points_line : public K_points
{
public:
	///\brief init followed by four steps read_k_points_coordinates(), read_reciprocal_basis(), calculate_k_points_coordinates() and output(). 
	///\return 0 if succeed. 
	int init();

	///\brief select bands at several k points. 
	///\param _selected_k_points the selected k points. 
	///\param _selected_num the number of k points selected. 
	///\return 0 if succeed. 
	int select_k_points(K_points &_selected_k_points, int _selected_num);

	///\brief select bands at several k points. 
	///\param _selected_k_points the selected k points. 
	///\param _index_of_k_points is a array of the index of k-points to be selected. 
	///\return 0 if succeed. 
	int select_k_points(K_points &_selected_k_points, Eigen::RowVectorXi &_index_of_k_points);

	K_points_line();///<\brief constructor. 

	~K_points_line();///<\brief destructor. 
private:
	Eigen::Matrix3d reciprocal_basis;///<\brief reciprocal basis. 
	Eigen::MatrixXd k_points_direct_coordinates;///<\brief select bands at several k points. 
	
	///\brief read k_points_direct_coordinates from ab-initio data. 
	///\return 0 if succeed. 
	int read_k_points_coordinates();
	
	///\brief read reciprocal basis from REPVEC. 
	///\return 0 if succeed. 
	int read_reciprocal_basis();
	
	///\brief calculate k points coordinates using k_points_direct_coordinates read from the ab-initio data. 
	///\return 0 if succeed. 
	int calculate_k_points_coordinates();
};

#endif  //_K_POINTS_LINE_H
