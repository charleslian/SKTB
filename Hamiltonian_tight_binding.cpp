#include "Hamiltonian_tight_binding.h"
#include "Atoms.h"
#include "K_points.h"
#include "Parameters.h"
#include <time.h>
#include <algorithm>

Hamiltonian_tight_binding::~Hamiltonian_tight_binding() 
{

}

Hamiltonian_tight_binding::Hamiltonian_tight_binding() 
{

}

int Hamiltonian_tight_binding::init(Atoms &_atoms, K_points &_k_points, Parameters &_parameters, Control &_control) 
{
	read_control(_control);//the max order of neighbors. 
	read_atoms(_atoms);//read the atoms. 
	read_k_points(_k_points);//read the k_points. 
	read_parameters(_parameters);//read the parameters. 
	calculate();
	//std::cout<<eigenvalue.eigenvalues;
	return 0;
}

int Hamiltonian_tight_binding::calculate()
{
	split_hamiltonian();
	int tight_binding_band_num=block_size*block_num;//read the number of the bands. 
	eigenvalue.eigenvalues=Eigen::MatrixXd(k_points_num,tight_binding_band_num);
	#pragma omp parallel for
	for (int k=0;k<k_points_num;k++)
	{
		//std::cout<<"calculate the "<<k<<"th points"<<std::endl;
		Eigen::MatrixXcd hamiltonian_temp(tight_binding_band_num,tight_binding_band_num);
		hamiltonian_temp=calculate_hamiltonian_of_k_index(k);
		calculate_eigenvalue(hamiltonian_temp,k);
	}
	return 0;
}

int Hamiltonian_tight_binding::split_hamiltonian()
{
	return 0;
}

int Hamiltonian_tight_binding::read_atoms(Atoms &_atoms)
{
	atoms=_atoms;
	return 0;
}

int Hamiltonian_tight_binding::read_control(Control & _control ) 
{
	control=_control;
	return 0;
}

int Hamiltonian_tight_binding::read_k_points(K_points &_k_points)
{
	k_points=_k_points;
	k_points_num=_k_points.k_points_num;
	return 0;
}

int Hamiltonian_tight_binding::read_parameters(Parameters &_parameters)
{
	parameters=_parameters;
	return 0;
}

Eigen::MatrixXcd Hamiltonian_tight_binding::calculate_self_energy(Atoms::Atom atom,int _k_points_index)
{
	Eigen::MatrixXcd hamiltonian_block(block_size,block_size);
	return hamiltonian_block;
}

Eigen::MatrixXcd Hamiltonian_tight_binding::calculate_hooping(Atoms::Atom atom1,Atoms::Atom atom2,int _k_points_index)
{
	Eigen::MatrixXcd hamiltonian_block(block_size,block_size);
	return hamiltonian_block;
}


Eigen::MatrixXcd Hamiltonian_tight_binding::calculate_hamiltonian_of_k_index(int _k_points_index) 
{
	Eigen::MatrixXcd hamiltonian_block(block_size,block_size);
	Eigen::MatrixXcd hamiltonian_of_k_index(block_size*block_num,block_size*block_num);
	hamiltonian_of_k_index.setZero();
	hamiltonian_block.setZero();
	for (int atom_index=0;atom_index<atoms.total_atom_num;atom_index++)
	{
		hamiltonian_block=calculate_self_energy(atoms.atom[atom_index],_k_points_index);
		hamiltonian_of_k_index.block(block_size*atom_index,block_size*atom_index,block_size,block_size)+=hamiltonian_block;
	}
	for (std::list<Atoms::Neighbor>::iterator neighbor_index = atoms.neighbors.begin(); neighbor_index != atoms.neighbors.end();neighbor_index++)
	{	
		Atoms::Atom atom1,atom2;
		int i=neighbor_index->index_i;
		int j=neighbor_index->index_j;
		atom1=atoms.atom[i];
		atom2=atoms.atom[j];
		atom1.coordinates=neighbor_index->coordinates1;
		atom2.coordinates=neighbor_index->coordinates2;
		hamiltonian_block=calculate_hooping(atom1,atom2,_k_points_index);
		hamiltonian_of_k_index.block(block_size*i,block_size*j,block_size,block_size)+=hamiltonian_block;
	}
	return hamiltonian_of_k_index;
}
int Hamiltonian_tight_binding::calculate_eigenvalue(Eigen::MatrixXcd hamiltonian,int _k_points_index)
{
	//sort the eigenvalues.
	int tight_binding_band_num=block_size*block_num;//read the number of the bands. 
	//solve the eigenfunction. 
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(hamiltonian);
	//sort the eigenvalues.
	double *tight_binding_eigen_temp=new double[tight_binding_band_num];
	for (int i=0;i<tight_binding_band_num;i++)
	{
		tight_binding_eigen_temp[i]=eigensolver.eigenvalues()(i);			
	}
	std::sort(tight_binding_eigen_temp,tight_binding_eigen_temp+tight_binding_band_num);
	for (int i=0;i<tight_binding_band_num;i++)
	{ 
		eigenvalue.eigenvalues(_k_points_index,i)=tight_binding_eigen_temp[i];
	}
	delete [] tight_binding_eigen_temp;
	return 0;
}

int Hamiltonian_tight_binding::calculate_DOS()
{
	//std::cout<<eigenvalue.eigenvalues.maxCoeff()<<std::endl;
	int N=control.number_of_point_in_dos_calculation;
	//std::cout<<N;
	Eigen::MatrixXd DOS(N,2);
	DOS.col(0).setLinSpaced(control.energy_interval[0],control.energy_interval[1]);
	for (int i=0;i<N;i++)
	{
		DOS(i,1)=fermi_smear(DOS(i,0),control.smear_sigma);
	}
	std::ofstream dos("DOS.dat");
	dos<<DOS;
	return 0;
}

double  Hamiltonian_tight_binding::fermi_smear(double energy,double ISMEAR)
{
	double dos=0;
	for (int i = 0; i < eigenvalue.eigenvalues.rows(); i++)
	{
		for (int j = 0; j <  eigenvalue.eigenvalues.cols(); j++)
		{
			double delta=fabs(energy-eigenvalue.eigenvalues(i,j))/ISMEAR;
			dos+=exp(-delta*delta)/(eigenvalue.eigenvalues.cols()*eigenvalue.eigenvalues.cols());
		}
	}
	return dos;
}
int Hamiltonian_tight_binding::output(char * filename)
{
	eigenvalue.output(filename);
	return 0;
}
