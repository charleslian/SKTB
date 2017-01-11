#include "Parameters_band.h"

Parameters_band::Parameters_band() 
{

}


Parameters_band::~Parameters_band() 
{

}

Parameters_band::Parameters_band(Atoms &_atoms) 
{
	int element_species=_atoms.element_species;
	int orbital_energy_species=2;//consider for self_energies:s, p.  
	int hooping_integral_species=4;//consider for hooping integrals:ss, sps, pps, ppp.  
	int parameter_num_per_hooping_integral=5;//each hooping integral has four parameters. 

	//read the orbital energies and hooping integrals.
	//for a system has N elements, X species of orbital-energy and Y species of hooping integrals.
	//there will be X*N orbital energies and Y*N*N hooping integrals in N elements system. 
	//each energy or hooping integral has parameter_num_per_hooping_integral parameters. 
	int num_of_lines=element_species*orbital_energy_species+element_species*element_species*hooping_integral_species;
	values=Eigen::MatrixXd(num_of_lines,parameter_num_per_hooping_integral);
	values.setZero();
	fitting_control=values;//the fitting_control has the same dimension of values.

	//read the values and fitting_control from file PARAM. 
	std::ifstream param("PARAM");
	//if the file PARAM not exists, set default values.
	if (!param)
	{
		double default_s_energy=-5.0;
		double default_p_energy=0.0;
		Eigen::Matrix4d hooping_parameters;
		hooping_parameters.setZero();
		hooping_parameters<<0,2.0,1.5,5.5,5.5,
			0,2.0,1.5,5.5,5.5,
			0,2.0,1.5,5.5,5.5,
			0,2.0,1.5,5.5,5.5;
		Eigen::MatrixXd delta(4,4);
		delta.setRandom();
		delta.col(0)*=10;
		hooping_parameters+=delta;
		int k=0;
		for (int i=0;i<element_species;i++)
		{
			values(k,0)=default_s_energy;
			values(k+1,0)=default_p_energy;
			k+=2;
		}
		for (int i=0;i<element_species*element_species;i++)
		{
			values.block(k,0,4,4)=hooping_parameters;
			k+=4;
		}
	}

	//if the file PARAM exists, read from it. 
	else
	{
		for (int i=0;i<num_of_lines;i++)
		{
			param>>values(i,0)>>values(i,1)>>values(i,2)>>values(i,3)>>values(i,4)>>fitting_control(i,0)>>fitting_control(i,1)>>fitting_control(i,2)>>fitting_control(i,3)>>fitting_control(i,4);
		}
	}
	param.close();

	//name the parameters. 
	parameter_name=new std::string[num_of_lines];
	int k=0;
	for (int i=0;i<element_species;i++)
	{
		parameter_name[k]="s orbital energy of element "+_atoms.element[i].name;
		parameter_name[k+1]="p orbital energy of element "+_atoms.element[i].name;
		k+=2;
	}
	for (int i=0;i<element_species;i++)
	{
		for (int j=0;j<element_species;j++)
		{
			parameter_name[k]="The ss between element "+_atoms.element[i].name+" and element "+_atoms.element[j].name;
			parameter_name[k+1]="The sps between element "+_atoms.element[i].name+" and element "+_atoms.element[j].name;
			parameter_name[k+2]="The pps between element "+_atoms.element[i].name+" and element "+_atoms.element[j].name;
			parameter_name[k+3]="The ppp between element "+_atoms.element[i].name+" and element "+_atoms.element[j].name;
			k+=4;
		}
	}

	//output to screen. 
	std::cout<<"---------------------------------------------------------------------------"<<std::endl;
	for (int i=0;i<num_of_lines;i++)
	{
		std::cout<<parameter_name[i]<<" is "<<values.row(i)<<std::endl;
	}

	std::ofstream param_vs_distance("DPARAM");
	Eigen::VectorXd distance;
	distance.setLinSpaced(100,0.1,10.1);
	Eigen::MatrixXd p(100,num_of_lines+1);
	p.col(0)=distance;
	for (int i=0;i<num_of_lines;i++)
	{
		double v0=values(i,0);
		double r0=values(i,1);
		double rc=values(i,2);
		double n=values(i,3);
		double nc=values(i,4);
		for (int j=0;j<100;j++)
		{
			p(j,i+1)=v0*pow(r0/distance(j),n)*exp(n*(-pow(distance(j)/rc,nc)+pow(r0/rc,nc)));
		}
	}
	param_vs_distance<<p;
}

