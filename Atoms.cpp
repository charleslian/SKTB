#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "Atoms.h"
#include <time.h>
Atoms::Atoms()
{

}

Atoms::~Atoms()
{

}

int Atoms::calculate_neighbor() 
{
	cutoff_distance_between_elements=Eigen::MatrixXd(element_species,element_species);
	cutoff_distance_between_elements.setOnes();
	cutoff_distance_between_elements*=6.0;

	std::ifstream read_cutoff("CUTOFF");
	for (int i=0;i<element_species;i++)
	{
		for (int j=0;j<element_species;j++)
		{
			read_cutoff>>cutoff_distance_between_elements(i,j);
			std::cout<<"The cutoff radius between element "<<element[i].name<<" with element "<<element[j].name<<" is "<<cutoff_distance_between_elements(i,j)<<std::endl;
		}
	}
	read_cutoff.close();
	
	//m n loop over all the atoms in the unit cell. 
	//output the neighbor distance to screen. 
	std::ofstream output_neighbor("NEIGH");
	Eigen::RowVector3d offset;
	offset.setZero();
	Atom atom_ghost;
	Eigen::RowVector3d delta_coordinates;
	double r;
	for(int m=0;m<total_atom_num;m++)
	{
		clock_t time1=clock();
		int element_m=atom[m].element.index;
		//std::list<Atom> neighbor_temp;
		for(int n=0;n<total_atom_num;n++)
		{
			int element_n=atom[n].element.index;
			int nth_neighbor=1+(int)(cutoff_distance_between_elements(element_m,element_n)/unit_cell_basis.row(0).norm());
			for(offset(0)=-nth_neighbor;offset(0)<=nth_neighbor;offset(0)++)
			{				
				for(offset(1)=-nth_neighbor;offset(1)<=nth_neighbor;offset(1)++)
				{
					for(offset(2)=-nth_neighbor;offset(2)<=nth_neighbor;offset(2)++)
					{
						//neighbor<<m<<" "<<n<<std::endl;
						//std::cout<<m<<" "<<n<<std::endl;
						//clock_t time1=clock();
						atom_ghost=atom[n];
						atom_ghost.coordinates.noalias()+=offset*unit_cell_basis;
						//calculate the distance
						delta_coordinates=atom_ghost.coordinates-atom[m].coordinates;
						r=delta_coordinates.norm();
						//if the distance is less than the cutoff radius, add this atom the the list of neighbor.  
						//clock_t time2=clock();
						if (r<cutoff_distance_between_elements(element_m,element_n)&&r>0.1)
						{
							Neighbor neighbor_temp;
							neighbor_temp.index_i=m;
							neighbor_temp.index_j=n;
							neighbor_temp.coordinates1=atom[m].coordinates;
							neighbor_temp.coordinates2=atom_ghost.coordinates;
							neighbors.push_back(neighbor_temp);
							//neighbor_temp.push_back(atom_ghost);
						}
						//clock_t time3=clock();
						//std::cout<<time1<<" "<<time2<<" "<<time3<<std::endl;
					}
				}
			}

		}
		clock_t time2=clock();
		//std::cout<<m<<" "<<time1<<" "<<time2<<" "<<std::endl;
		//output the neighbor list to screen
		//output_neighbor<<"The neighbor coordinates of atom "<<m<<std::endl;
	}		
	for (std::list<Neighbor>::iterator i = neighbors.begin(); i != neighbors.end(); ++i)
	{
		//std::cout <<i->index_i<<" "<<i->index_j<<" "<<i->distance<<std::endl;
		output_neighbor<<i->index_i<<" "<<i->index_j<<std::endl<<i->coordinates1<<std::endl<<i->coordinates2<<std::endl;
	}
	output_neighbor<<std::endl;
	output_neighbor.close();
	return 0;
}

int Atoms::calculate_reciprocal_basis() 
{
	const double PI=3.1415926;
	double unit_cell_volume=unit_cell_basis.row(0).dot(unit_cell_basis.row(1).cross(unit_cell_basis.row(2)));

	//calculate the reciprocal basis. 
	reciprocal_basis.row(0)=2*PI/unit_cell_volume*unit_cell_basis.row(1).cross(unit_cell_basis.row(2));
	reciprocal_basis.row(1)=2*PI/unit_cell_volume*unit_cell_basis.row(2).cross(unit_cell_basis.row(0));
	reciprocal_basis.row(2)=2*PI/unit_cell_volume*unit_cell_basis.row(0).cross(unit_cell_basis.row(1));
	std::cout<<"Reciprocal basis:"<<std::endl<<reciprocal_basis<<std::endl;

	//output the reciprocal basis to REPVEC
	std::ofstream repvec("REPVEC");
	repvec<<reciprocal_basis;
	repvec.close();

	return 0;
}