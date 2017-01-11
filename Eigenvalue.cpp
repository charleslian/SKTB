#include <fstream>
#include <iostream>
#include <algorithm>
#include "Eigen/Eigenvalues"
#include "Eigenvalue.h"

Eigenvalue::Eigenvalue()
{

}

Eigenvalue::~Eigenvalue()
{

}


int Eigenvalue::select_bands(Eigenvalue &_selected_bands, int &_lowest_band, int &_highest_band) 
{
	Eigen::MatrixXd eigenvalues_temp(eigenvalues.rows(),_highest_band-_lowest_band);
	for(int j=_lowest_band;j<_highest_band;j++)
	{
		eigenvalues_temp.col(j-_lowest_band)=eigenvalues.col(j);
	}
	_selected_bands.eigenvalues=eigenvalues_temp;
	return 0;
}

int Eigenvalue::select_bands(Eigenvalue &_selected_bands, Eigen::RowVectorXi &_index_of_bands)
{
	Eigen::MatrixXd eigenvalues_temp(eigenvalues.rows(),_index_of_bands.cols());
	for(int j=0;j<eigenvalues_temp.cols();j++)
	{
		eigenvalues_temp.col(j)=eigenvalues.col(_index_of_bands(j));
	}
	_selected_bands.eigenvalues=eigenvalues_temp;
	return 0;
}

int Eigenvalue::select_k_points(Eigenvalue &_selected_k_points, int &_selected_num) 
{
	int kspace=eigenvalues.rows()/_selected_num;
	_selected_k_points.eigenvalues=Eigen::MatrixXd(_selected_num,eigenvalues.cols());
	for (int i=0;i<_selected_num;i++)
	{
		_selected_k_points.eigenvalues.row(i)=eigenvalues.row(i*kspace);
	}
	return 0;
}

int Eigenvalue::select_k_points(Eigenvalue &_selected_k_points, Eigen::RowVectorXi &_index_of_k_points) 
{
	_selected_k_points.eigenvalues=Eigen::MatrixXd(_index_of_k_points.cols(),eigenvalues.cols());
	for (int i=0;i<_selected_k_points.eigenvalues.cols();i++)
	{
		_selected_k_points.eigenvalues.row(i)=eigenvalues.row(_index_of_k_points(i));
	}
	return 0;
}

int Eigenvalue::output(char *filename) 
{	
	band_x_coordinates=Eigen::VectorXd (eigenvalues.rows());
	band_x_coordinates.setZero();
	//read band_x_coordinates from k_coordinates.dat. 
	std::ifstream k_coordinates("k_coordinates.dat");
	for (int i=0;i<band_x_coordinates.rows();i++)
	{
		double temp;
		k_coordinates>>band_x_coordinates(i)>>temp>>temp>>temp;
	}
	k_coordinates.close();

	//set the first column of eigenvalues to be band_x_coordinates. 
	Eigen::MatrixXd eigenvalue_out(eigenvalues.rows(),eigenvalues.cols()+1);
	eigenvalue_out.col(0)=band_x_coordinates;
	//std::cout<<eigenvalue_out.cols();
	for(int i=1;i<eigenvalue_out.cols();i++)
	{
		eigenvalue_out.col(i)=eigenvalues.col(i-1);
	}

	//output data to filename. 
	std::ofstream energy_band(filename);
	energy_band<<"k Energy"<<" ";
	for (int i=2;i<eigenvalue_out.cols();i++)
	{
		energy_band<<i<<" ";
	}
	energy_band<<std::endl;
	energy_band<<"A eV"<<std::endl;
	energy_band<<eigenvalue_out;
	energy_band.close();

	return 0;
}

int Eigenvalue::calculate_band_gap_and_Vf(int _electron_number,double smear)
{
	int highest_occupied_band=_electron_number/2-1;
	//std::cout<<highest_occupied_band<<std::endl;
	//std::cout<<band.col(highest_occupied_band)<<std::endl;
	Eigen::VectorXd HOMO=eigenvalues.col(highest_occupied_band);
	Eigen::VectorXd LUMO=eigenvalues.col(highest_occupied_band+1);
	// 	std::cout<<HOMO<<std::endl;
	// 	std::cout<<LUMO<<std::endl;
	// 	std::ofstream test_out("Dirac.dat");
	// 	for (int i = 0; i < x_coordinates.rows(); i++)
	// 	{
	// 		test_out<<" "<<x_coordinates(i)<<" "<<HOMO(i)<<" "<<LUMO(i)<<std::endl;
	// 	}
	//std::cout<<" "<<x_coordinates(i)<<" "<<HOMO(i)<<" "<<LUMO(i)<<std::endl;
	unsigned int bottom_index,top_index;
	double gap_bottom=HOMO.maxCoeff(&bottom_index);
	double gap_top=LUMO.minCoeff(&top_index);
	
	
	Eigen::MatrixXd drift(eigenvalues.rows(),eigenvalues.cols());
	drift.setOnes();
	drift.col(0).setZero();
	eigenvalues-=drift*gap_bottom;
	this->output("tight-binding-minus-Ef.dat");

	std::cout<<"Band gap is "<<gap_top-gap_bottom<<std::endl;
	std::cout<<"The top of band gap is at "<<top_index<<", the bottom is at "<<bottom_index<<std::endl;
	std::ofstream output_band_gap("band_gap.dat");
	output_band_gap<<gap_top-gap_bottom<<std::endl;
	output_band_gap.close();

	int interval=0;
	double energy_difference=0;
	double average_velocity=0;
	for (int i = 1; i < eigenvalues.rows()-1; i++)
	{
		energy_difference=HOMO(i)-HOMO(bottom_index);
		//std::cout<<"test:"<<energy_difference<<std::endl;
		if (fabs(energy_difference)<smear&&fabs(energy_difference)>0)
		{
			double length1=band_x_coordinates(i)-band_x_coordinates(i+1);
			double length2=band_x_coordinates(i)-band_x_coordinates(i-1);
			double velocity=0;
			if (fabs(length1)>0.001)
			{
				velocity+=fabs((HOMO(i)-HOMO(i+1))/(length1));
				interval++;
			}
			if (fabs(length2)>0.001)
			{
				velocity+=fabs((HOMO(i)-HOMO(i-1))/(length2));
				interval++;
			}
			average_velocity+=velocity;
		}
	}
	//std::cout<<interval<<std::endl;
	if (interval==0)
	{
		average_velocity=0;
	}
	else
	{
		average_velocity/=interval;
	}
	std::ofstream output_velocity("Fermi_velocity.dat");
	average_velocity*=(1.6/1.054);
	std::cout<<"The Fermi velocity is "<<average_velocity<<std::endl;
	output_velocity<<average_velocity<<std::endl;
	//std::cout<<band<<std::endl;
	return 0;
}
