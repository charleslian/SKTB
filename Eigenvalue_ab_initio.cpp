#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "Eigen/Eigenvalues"
#include "Eigenvalue_ab_initio.h"

Eigenvalue_ab_initio::Eigenvalue_ab_initio() 
{

}

Eigenvalue_ab_initio::~Eigenvalue_ab_initio() 
{

}

int Eigenvalue_ab_initio::calculate_eigenvalue()
{
	int k_points_num;//the eigenvalues calculated. 
	int band_num;//the eigenvalues calculated. 
	int ispin;//ispin=0:not polarized, ispin=1:polarized. 
	double fermi_energy;//the eigenvalues calculated. 
	std::string temp;
	double temp1;
	
	std::ifstream doscar("DOSCAR");//read Fermi energy from DOSCAR file. 
	//discard useless data. 
	for (int i=0;i<5;i++)
	{
		std::getline(doscar,temp);
	}
	doscar>>temp1>>temp1>>temp1>>fermi_energy>>temp1;//read Fermi energy.
	doscar.close();
	
	std::ifstream eigenval("EIGENVAL");//read eigenvalues from EIGENVAL file. 
	eigenval>>temp1>>temp1>>temp1>>ispin;//read ispin. 
	//discard useless data. 
	for (int i=0;i<5;i++)
	{
		std::getline(eigenval,temp);
	}
	eigenval>>temp1>>k_points_num>>band_num;//read total number of k-points and the number of bands.
	eigenvalues=Eigen::MatrixXd(k_points_num,band_num);//read eigenvalues. 
	for (int i=0;i<k_points_num;i++)
	{
		//read three blank lines. 
		for (int k=0;k<3;k++)
		{
			std::getline(eigenval,temp);
		}
		for (int j=0;j<band_num;j++)
		{
			if (ispin==2)
			{
				eigenval>>temp1>>eigenvalues(i,j)>>temp1;
			}
			else
			{
				eigenval>>temp1>>eigenvalues(i,j);
			}
			eigenvalues(i,j)-=fermi_energy;
		}
	}
	eigenval.close();

	//output to screen. 
	std::cout<<"----------------------ab initio data---------------------------------------"<<std::endl;
	std::cout<<"Fermi energy is:"<<fermi_energy<<std::endl;
	if (ispin==2)
	{
		std::cout<<"Calculation is spin-polarized."<<std::endl;
	} 
	else
	{
		std::cout<<"Calculation is not spin-polarized."<<std::endl;
	}
	std::cout<<"The number of calculated band is "<<band_num<<std::endl;
	std::cout<<"The number of calculated k points is "<<k_points_num<<std::endl;
	return 0;
}


