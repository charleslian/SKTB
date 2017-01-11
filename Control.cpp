#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "Control.h"

Control::Control() 
{

}

Control::Control(char * filename) 
{
	//set default value
	set_default_value();

	//read the control file
	read_control_tag(filename);

	//output the information. 
	output();
}

Control::~Control() 
{

}

int Control::set_default_value()
{
	fitting_needed=false;
	SOC_calculation_needed=false;
	split_needed=false;
	total_energy_calculation_needed=false;
	band_calculation_needed=false;
	total_energy_calculation_needed=false;
	lanczos_calculation_needed=false;
	DOS_calculation_needed=false;

	order_of_lanzos_matrix=10;
	number_of_point_in_dos_calculation=100;
	eta=0.001;
	energy_interval[0]=-2;energy_interval[1]=2;

	k_mesh_num(0)=k_mesh_num(1)=k_mesh_num(2)=1;

	criteria=0.00001;
	num_of_fitting_point=1;
	max_neighbor_order_for_band_calculation=1;
	k_mesh_num(0)=1,k_mesh_num(1)=1,k_mesh_num(2)=1;
	fitting_band[0]=1;fitting_band[1]=1;
	electric_field_direction=2;
	electric_field_strength=0;
	return 0;
};

int Control::read_control_tag(char * filename)
{
	std::ifstream control(filename);
	std::string temp;
	while (!control.eof())
	{
		std::getline(control,temp);

		//find the split_needed input. 
		if (temp.find("ISPL")!=temp.npos)
		{
			if (temp.find("1")!=temp.npos)
			{
				split_needed=true;
			}
			else
			{
				split_needed=false;
			}
		}

		if (temp.find("IDOS")!=temp.npos)
		{
			if (temp.find("1")!=temp.npos)
			{
				DOS_calculation_needed=true;
			}
			else
			{
				DOS_calculation_needed=false;
			}
		}
		//find the fitting_needed input. 
		if (temp.find("IFIT")!=temp.npos)
		{
			if (temp.find("1")!=temp.npos)
			{
				fitting_needed=true;
			}
			else
			{
				fitting_needed=false;
			}
		}

		//find the SOC_calculation_needed input. 
		if (temp.find("ISOC")!=temp.npos)
		{
			if (temp.find("1")!=temp.npos)
			{
				SOC_calculation_needed=true;
			}
			else
			{
				SOC_calculation_needed=false;
			}
		}

		//find the total_energy_calculation_needed input. 
		if (temp.find("IBAN")!=temp.npos)
		{
			if (temp.find("1")!=temp.npos)
			{
				band_calculation_needed=true;
			}
			else
			{
				band_calculation_needed=false;
			}
		}

		//find the criteria input. 
		if (temp.find("DCRI")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>criteria;
		}

		//find the num_of_fitting_point input. 
		if (temp.find("INKP")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>num_of_fitting_point;
		}

		if (temp.find("ISMEAR")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>smear_sigma;
		}

		if (temp.find("DDIP")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>electric_field_direction;
		}

		if (temp.find("EDIP")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>electric_field_strength;
		}

		if (temp.find("IENC")!=temp.npos)
		{
			if (temp.find("1")!=temp.npos)
			{
				total_energy_calculation_needed=true;
			}
			else
			{
				total_energy_calculation_needed=false;
			}
		}
		//find the max input. 
		if (temp.find("INNO")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>max_neighbor_order_for_band_calculation;
		}

		//find the fitting_band input. 
		if (temp.find("IMMB")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>fitting_band[0]>>fitting_band[1];
		}

		//find the k_mesh_num input. 
		if (temp.find("IKMG")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			temp.substr(position,position);
			iss>>k_mesh_num(0)>>k_mesh_num(1)>>k_mesh_num(2);
		}

		//find the num_of_fitting_point input. 
		if (temp.find("OOLM")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>order_of_lanzos_matrix;
		}

		//find the max input. 
		if (temp.find("NDOS")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>number_of_point_in_dos_calculation;
		}

		if (temp.find("ETA")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>eta;
		}

		//find the fitting_band input. 
		if (temp.find("EMTM")!=temp.npos)
		{
			std::string::size_type position;
			position=temp.find("=");
			position++;
			std::istringstream iss(temp.substr(position));
			iss>>energy_interval[0]>>energy_interval[1];
		}
	}
	return 0;
}
int Control::output()
{
	
	std::cout<<"----------------------control tag------------------------------------------"<<std::endl;
	
	if (fitting_needed)
	{
		std::cout<<"The parameters will be adjusted to fit the ab initio data."<<std::endl;
	} 
	else
	{
		std::cout<<"The parameters will not be changed."<<std::endl;
	}

	if (band_calculation_needed)
	{
		std::cout<<"The band will be calculated."<<std::endl;
	} 
	else
	{
		std::cout<<"The band will not be calculated."<<std::endl;
	}

	if (total_energy_calculation_needed)
	{
		std::cout<<"The total energy will be calculated."<<std::endl;
	} 
	else
	{
		std::cout<<"The total energy will not be calculated."<<std::endl;
	}

	std::cout<<"The fitting criteria is "<<criteria<<std::endl;
	std::cout<<"The number of k points used for fitting is "<<num_of_fitting_point<<std::endl;
	std::cout<<"The number of bands used for fitting is from "<<fitting_band[0]<<" to "<<fitting_band[1]<<std::endl;
	std::cout<<"The range of repeated cell is "<<-max_neighbor_order_for_band_calculation<<" to "<<max_neighbor_order_for_band_calculation<<std::endl;
	std::cout<<"The K mesh grid is "<<k_mesh_num(0)<<" x "<<k_mesh_num(1)<<" x "<<k_mesh_num(2)<<std::endl;
	
	return 0;
}