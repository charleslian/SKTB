#include "K_points_grid.h"

K_points_grid::K_points_grid()
{

}

K_points_grid::~K_points_grid()
{

}

K_points_grid::K_points_grid(Eigen::RowVector3d &_k_mesh_num)
{
	read_reciprocal_basis();
	calculate_k_grid(_k_mesh_num);
	output();
}


int K_points_grid::read_reciprocal_basis()
{
	//read reciprocal basis from REPVEC file. 
	std::ifstream repvec("REPVEC");
	for (int i=0;i<3;i++)
	{
		repvec>>reciprocal_basis(i,0)>>reciprocal_basis(i,1)>>reciprocal_basis(i,2);
	}
	repvec.close();
	return 0;
}

int K_points_grid::calculate_k_grid(Eigen::RowVector3d &_k_mesh_num)
{
	//the number of k points is the number of grid knots. 
	k_points_num=(int)(_k_mesh_num(0)*_k_mesh_num(1)*_k_mesh_num(2));
	
	//calculate the k-points coordinates. 
	k_points_coordinates=Eigen::MatrixXd(k_points_num,3);
	k_points_coordinates.setZero();
	Eigen::RowVector3i k_index;
	k_index.setZero();
	int index=0;
	for (k_index(0)=0;k_index(0)<_k_mesh_num(0);k_index(0)++)
	{
		for (k_index(1)=0;k_index(1)<_k_mesh_num(1);k_index(1)++)
		{
			for (k_index(2)=0;k_index(2)<_k_mesh_num(2);k_index(2)++)
			{		
				for (int i=0;i<3;i++)
				{	
					k_points_coordinates.row(index)+=(2*k_index(i)-_k_mesh_num(i)+1)/(2*_k_mesh_num(i))*reciprocal_basis.row(i);
				}
				index++;
			}
		}
	}

	if (0)
	{
		Eigen::Matrix3d rotation;rotation.setZero();
		double theta_in_angle;
		theta_in_angle=180;
		double theta=theta_in_angle/180.0*3.1415926;
		std::cout<<sin(theta)<<" "<<cos(theta)<<std::endl;
		rotation(0,0)=rotation(1,1)=cos(theta);
		rotation(0,1)=sin(theta);rotation(1,0)=-sin(theta);
		rotation(2,2)=1;

		Eigen::MatrixXd k_points_coordinates_temp=k_points_coordinates*rotation;
		Eigen::VectorXi exist(k_points_num);exist.setZero();
		for (int i=0;i<k_points_num;i++)
		{
			for (int j = 0; j < k_points_num; j++)
			{
				if (abs((k_points_coordinates_temp.row(i)-k_points_coordinates.row(j)).norm())<0.01)
				{
					std::cout<<i<<" "<<j<<std::endl;
					exist(i)=1;
					break;
				}
			}
		}
		std::cout<<exist<<std::endl;
		if (exist.minCoeff()==1)
		{
			std::cout<<"The rotation "<<theta_in_angle<<" is conserved;"<<std::endl;
		}
		//std::cout<<(k_points_coordinates_temp-k_points_coordinates).norm()<<std::endl;

	}

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
