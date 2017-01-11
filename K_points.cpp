#include "K_points.h"

K_points::K_points() 
{

}

K_points::~K_points() 
{

}

int K_points::output() 
{
	//output the k-points coordinates to file k_coordinates.dat. 
	std::ofstream k_coordinates("k_coordinates.dat");
	for (int i=0;i<k_points_coordinates.rows();i++)
	{
		k_coordinates<<band_x_coordinates(i)<<" ";
		for (int j=0;j<k_points_coordinates.cols();j++)
		{
			k_coordinates<<k_points_coordinates(i,j)<<" ";
		}
		k_coordinates<<std::endl;
	}
	k_coordinates.close();
	return 0;
}



