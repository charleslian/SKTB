///\brief Define Eigenvalue class. 
///\details class to store eigenvalues and select eigenvalues at certain bands or k-points. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved. 


#if !defined(_EIGENVALUE_H)
#define _EIGENVALUE_H

#include "Eigen/Dense"

class Eigenvalue
{
public:
	Eigen::MatrixXd eigenvalues;///<\brief the eigenvalues calculated. 
	Eigen::VectorXd band_x_coordinates;
	Eigen::MatrixXcd *eigenvectors;
	int calculate_band_gap_and_Vf(int _electron_number,double smear);
	///\brief selected some bands from eigenvalues. 
	///\param _selected_bands is another Eigenvalue type data store the selected eigenvalues.  
	///\param _lowest_band is the lowest band to be selected. 
	///\param _highest_band is the highest band to be selected. 
	///\return 0 if succeed. 
	int select_bands(Eigenvalue &_selected_bands, int &_lowest_band, int &_highest_band);

	///\brief selected some bands from eigenvalues. 
	///\param _selected_bands is another Eigenvalue type data store the selected eigenvalues.  
	///\param _index_of_bands is a array of the index of bands to be selected.  
	///\return 0 if succeed. 
	int select_bands(Eigenvalue &_selected_bands, Eigen::RowVectorXi &_index_of_bands);

	///\brief selected some k-points from eigenvalues, the k-points is selected uniquely. 
	///\param _selected_k_points is another Eigenvalue type data store the selected eigenvalues.  
	///\param _selected_num is number of k-points to be selected. 
	///\return 0 if succeed. 
	int select_k_points(Eigenvalue &_selected_k_points, int &_selected_num);

	///\brief selected some k-points from eigenvalues, the k-points to be selected is manually specified.  
	///\param _selected_k_points is another Eigenvalue type data store the selected eigenvalues.  
	///\param _index_of_k_points is a array of the index of k-points to be selected.  
	///\return 0 if succeed. 
	int select_k_points(Eigenvalue &_selected_k_points, Eigen::RowVectorXi &_index_of_k_points);

	///\brief output the eigenvalues. 
	///\param *filename is the output filename. 
	///\return 0 if succeed. 
	int output(char *filename);

	
	Eigenvalue();///<\brief constructor. 

	~Eigenvalue();///<\brief destructor. 

};

#endif  //_EIGENVALUE_H
