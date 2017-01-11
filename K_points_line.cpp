#include "K_points_line.h"
K_points_line::K_points_line()
{

}

K_points_line::~K_points_line()
{

}

int K_points_line::init()
{
	read_k_points_coordinates();//read k points direct coordinates from file EIGENVAL. 
	read_reciprocal_basis();//read reciprocal basis file REPVEC. 
	calculate_k_points_coordinates();//calculate k points coordinates using direct coordinates. 
	output();

	return 0;
}

int K_points_line::select_k_points(K_points &_selected_k_points, int _selected_num)
{
	Eigen::MatrixXd k_points_coordinates_temp;
	int kspace=k_points_num/_selected_num;
	_selected_k_points.k_points_num=_selected_num;
	_selected_k_points.k_points_coordinates=Eigen::MatrixXd(_selected_num,3);
	for (int i=0;i<_selected_num;i++)
	{
		_selected_k_points.k_points_coordinates.row(i)=k_points_coordinates.row(i*kspace);
	}
	return 0;
}

int K_points_line::select_k_points(K_points &_selected_k_points, Eigen::RowVectorXi &_index_of_k_points)
{
	Eigen::MatrixXd k_points_coordinates_temp;
	_selected_k_points.k_points_num=_index_of_k_points.cols();
	_selected_k_points.k_points_coordinates=Eigen::MatrixXd(_index_of_k_points.cols(),3);
	for (int i=0;i<_index_of_k_points.cols();i++)
	{
		_selected_k_points.k_points_coordinates.row(i)=k_points_coordinates.row(_index_of_k_points(i));
	}
	return 0;
}

int K_points_line::read_k_points_coordinates()
{
	int band_num;//the eigenvalues calculated. 
	std::string temp;
	double temp1;
	std::ifstream eigenval("EIGENVAL");//read eigenvalues from EIGENVAL file. 
	//discard useless data. 
	for (int i=0;i<5;i++)
	{
		std::getline(eigenval,temp);
	}
	eigenval>>temp1>>k_points_num>>band_num;//read total number of k-points and the number of bands.
	k_points_direct_coordinates=Eigen::MatrixXd(k_points_num,3);//read eigenvalues. 
	for (int i=0;i<k_points_num;i++)
	{
		//read three blank lines. 
		for (int k=0;k<2;k++)
		{
			std::getline(eigenval,temp);
		}
		eigenval>>k_points_direct_coordinates(i,0)>>k_points_direct_coordinates(i,1)>>k_points_direct_coordinates(i,2)>>temp1;
		for (int j=0;j<band_num;j++)
		{
			std::getline(eigenval,temp);
		}
	}
	eigenval.close();
	return 0;
}

int K_points_line::read_reciprocal_basis()
{
	std::ifstream repvec("REPVEC");
	for (int i=0;i<3;i++)
	{
		repvec>>reciprocal_basis(i,0)>>reciprocal_basis(i,1)>>reciprocal_basis(i,2);
	}
	repvec.close();
	return 0;
}

int K_points_line::calculate_k_points_coordinates()
{
	//calculate the k-points coordinates. 
	k_points_coordinates=Eigen::MatrixXd(k_points_num,3);
	k_points_coordinates=k_points_direct_coordinates*reciprocal_basis;
	
	//calculate the x coordinates of energy band. 
	band_x_coordinates=Eigen::RowVectorXd(k_points_num);
	band_x_coordinates(0)=0.0;
	for (int i=1;i<k_points_num;i++)
	{
		Eigen::RowVector3d delta_k_points_coordinates=k_points_coordinates.row(i)-k_points_coordinates.row(i-1);
		band_x_coordinates(i)=band_x_coordinates(i-1)+delta_k_points_coordinates.norm();
	}
	return 0;
}

