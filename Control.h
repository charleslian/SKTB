///\brief Define Control class. 
///\details class used by the main function to control the process. 
///\date  2013/12/10. 
///\author  Chao Lian. 
///\copyright Copyright 2013 Chao Lian. All rights reserved.

#if !defined(_CONTROL_H)
#define _CONTROL_H

#include "Eigen/Dense"

class Control
{
public:

	int order_of_lanzos_matrix;///<\brief the number of the k-points used for fitting. 
	double energy_interval[2];///<\brief the indices of the lowest band and the highest band used for fitting. 
	double eta;///<\brief the order of the neighbors used for mobius calculation.  
	int number_of_point_in_dos_calculation;///<\brief the order of the neighbors used for band calculation. 
	double smear_sigma;

	bool total_energy_calculation_needed;///<\brief whether the total-energy calculation is needed. 
	bool fitting_needed;///<\brief whether the fitting of parameters is needed. 
	bool split_needed;///<\brief whether the split of p-orbital is considered. 
	bool band_calculation_needed;///<\brief whether the band-calculation is needed. 
	bool SOC_calculation_needed;///<\brief whether the SOC_band-calculation is needed. 
	bool lanczos_calculation_needed;///<\brief whether the SOC_band-calculation is needed. 
	bool DOS_calculation_needed;//
	int num_of_fitting_point;///<\brief the number of the k-points used for fitting. 
	int fitting_band[2];///<\brief the indices of the lowest band and the highest band used for fitting. 
	double max_neighbor_order_for_band_calculation;///<\brief the order of the neighbors used for band calculation. 
	double criteria;///<\brief the fitting criteria. 
	double electric_field_strength;
	int electric_field_direction;
	Eigen::RowVector3d k_mesh_num;///<\brief the number of mesh-grid points in x y z direction. 

	Control();///<\brief constructor. 

	~Control();///<\brief destructor. 
	Control(char * filename) ;
	int set_default_value();
	int read_control_tag(char * filename);
	int output();
};

#endif  //_CONTROL_H