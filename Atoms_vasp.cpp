#include <fstream>
#include <iostream>
#include <string>
#include <cctype>
#include "Atoms_vasp.h"

Atoms_vasp::Atoms_vasp() 
{
	read();
	calculate_neighbor();
	calculate_reciprocal_basis();
}

Atoms_vasp::~Atoms_vasp() 
{
	delete [] atom;
	delete [] element;
}

int Atoms_vasp::read() 
{
	double universal_scaling_factor;//lattice constant. 
	std::string system_name;//system name. 
	std::string	unit_cell_basis_type;//basis type: Cartesian or Direct.
	std::string	temp;//temp variable to read useless data. 
	
	//read from POSCAR file. 
	std::ifstream poscar("POSCAR");
	
	//read the system name and lattice constant. 
	std::getline(poscar,system_name);
	poscar>>universal_scaling_factor;

	//read the basis. 
	for(int i=0;i<3;i++)
	{
		poscar>>unit_cell_basis(i,0)>>unit_cell_basis(i,1)>>unit_cell_basis(i,2);
	}
	unit_cell_basis*=universal_scaling_factor;
	
	//read the names of elements. 
	std::getline(poscar,temp);
	std::getline(poscar,temp);

	//calculate the total element species. 
	element_species=0;
	for (unsigned int i=0;i<temp.length();i++)
	{
		if(std::isalpha(temp[i]))
		{
			if (i==0||(!std::isalpha(temp[i-1])))
			{
				element_species++;
			}
		}
	}

	//read the names of elements. 
	element=new Element[element_species];
	int	j=0;
	for (unsigned int i=0;i<temp.length();i++)
	{
		if(std::isalpha(temp[i]))
		{
			if (i==0||(!std::isalpha(temp[i-1])))
			{
				element[j].name=temp[i];
				j++;
			}
			else
			{
				element[j-1].name+=temp[i];
			}
		}
	}

	//read the number of atoms per element. 
	Eigen::RowVectorXi atom_num_of_element(element_species);
	for(int i=0;i<element_species;i++)
	{ 
		poscar>>atom_num_of_element(i);
	}

	//calculate the total number of atoms. 
	total_atom_num=0;
	for(int i=0;i<element_species;i++) 
	{
		total_atom_num+=atom_num_of_element(i);
	}
	
	//read the type of basis. 
	std::getline(poscar,temp);
	std::getline(poscar,unit_cell_basis_type);

	//treat Cartesian and Direct type separately. 
	//read the coordinates and convert to Cartesian coordinates. 
	atom=new Atoms::Atom[total_atom_num];
	Eigen::MatrixXd direct_coordinates(total_atom_num,3);
	if(unit_cell_basis_type[0]=='D'||unit_cell_basis_type[0]=='d')
	{
		for(int i=0;i<total_atom_num;i++)
		{
			poscar>>direct_coordinates(i,0)>>direct_coordinates(i,1)>>direct_coordinates(i,2);//
			atom[i].coordinates=direct_coordinates.row(i)*unit_cell_basis;//
			atom[i].index=i;//
		}
	}
	else
	{
		for(int i=0;i<total_atom_num;i++)
		{
			poscar>>atom[i].coordinates(0)>>atom[i].coordinates(1)>>atom[i].coordinates(2);
			atom[i].index=i;
		}
	}

	//set the element of every atom
	int k=0;	
	for(int i=0;i<element_species;i++)
	{
		for(int j=0;j<atom_num_of_element(i);j++)
		{
			atom[k].element.index=i;
			atom[k].element.name=element[i].name;
			k++;
		}	
	}


	//output the information. 
	std::cout<<"---------------------------------------------------------------------------"<<std::endl;
	std::cout<<"System name: "<<system_name<<std::endl;
	std::cout<<"lattice constant: "<<universal_scaling_factor<<std::endl;
	std::cout<<"---------------------------------------------------------------------------"<<std::endl;
	std::cout<<"Unit cell basis:"<<std::endl<<unit_cell_basis<<std::endl;
	std::cout<<"atom numbers of elements:"<<atom_num_of_element<<std::endl;
	std::cout<<"total number of atoms:"<<total_atom_num<<std::endl;
	std::cout<<"The basis type is "<<unit_cell_basis_type<<std::endl;
	std::cout<<"---------------------------------------------------------------------------"<<std::endl;
	for(int i=0;i<element_species;i++) 
	{
		std::cout<<"Element "<<element[i].name<<" atom number is "<<atom_num_of_element(i)<<std::endl;
	}

	std::cout<<"---------------------------------------------------------------------------"<<std::endl;
	for(int i=0;i<total_atom_num;i++)
	{
		std::cout<<"Atom:"<<i<<" ";
		std::cout<<"Element:"<<atom[i].element.name<<" ";
		std::cout<<"coordinates:"<<atom[i].coordinates<<std::endl;
	}
	poscar.close();

	total_electron_num=0;
	for (int i=0;i<element_species;i++)
	{
		if (element[i].name=="Si")
		{
			element[i].electron_number=4;
		}
		else if (element[i].name=="C")
		{
			element[i].electron_number=4;
		}
		else if (element[i].name=="B")
		{
			element[i].electron_number=3;
		}
		else if (element[i].name=="Al")
		{
			element[i].electron_number=3;
		}
		else if (element[i].name=="Tl")
		{
			element[i].electron_number=3;
		}
		else if (element[i].name=="P")
		{
			element[i].electron_number=5;
		}
		total_electron_num+=element[i].electron_number*atom_num_of_element(i);
	}
	std::cout<<"total_electron_num: "<<total_electron_num<<std::endl;
	
	return 0;
}

