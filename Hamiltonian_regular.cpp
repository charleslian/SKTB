#include "Hamiltonian_regular.h"
#include "Atoms.h"
#include "K_points.h"
#include "Parameters.h"

Hamiltonian_regular::Hamiltonian_regular()
{

}

Hamiltonian_regular::~Hamiltonian_regular()
{

}

Hamiltonian_regular::Hamiltonian_regular(Atoms &_atoms, K_points &_k_points, Parameters &_parameters, Control &_control) 
{
	init(_atoms, _k_points, _parameters,_control);
}

int Hamiltonian_regular::split_hamiltonian() 
{
	//four orbitals for one block. 
	block_size=4;

	//block number is equal to the number of total atoms.
	block_num=atoms.total_atom_num;

	return 0;
}

 
Eigen::MatrixXcd Hamiltonian_regular::calculate_self_energy(Atoms::Atom atom,int _k_points_index)
{
	//read the orbital energy. 
	int orbital_energy_species=2;//orbital energies:s, p. 
	Eigen::MatrixXd orbital_energy(orbital_energy_species,parameters.values.cols());
	for (int i=0;i<orbital_energy_species;i++)
	{
		orbital_energy.row(i)=parameters.values.row(i+orbital_energy_species*atom.element.index);
	}

	//calculate the element of self_energy. 
	Eigen::MatrixXcd hamiltonian_block(block_size,block_size);
	hamiltonian_block.setZero();
	hamiltonian_block(0,0)=orbital_energy(0,0)+atom.coordinates(control.electric_field_direction)*control.electric_field_strength;
	hamiltonian_block(1,1)=orbital_energy(1,0)+atom.coordinates(control.electric_field_direction)*control.electric_field_strength;
	hamiltonian_block(2,2)=orbital_energy(1,0)+atom.coordinates(control.electric_field_direction)*control.electric_field_strength;
	hamiltonian_block(3,3)=orbital_energy(1,0)+atom.coordinates(control.electric_field_direction)*control.electric_field_strength;
	//std::cout<<hamiltonian_block<<std::endl;
	return hamiltonian_block;
}

Eigen::MatrixXcd Hamiltonian_regular::calculate_hooping(Atoms::Atom atom1,Atoms::Atom atom2,int _k_points_index)
{
	//calculate the displacement and azimuth cosine. 
	Eigen::RowVector3d delta_coordinates=atom2.coordinates-atom1.coordinates;
	double r=delta_coordinates.norm();
	Eigen::RowVector3d Azimuth_cos=delta_coordinates/r;

	//read hooping integrals parameters and calculate the hooping integrals. 
	int hooping_integral_species=4;//hooping integrals: ss, sp, pps, ppp. 
	Eigen::MatrixXd hooping_integral(hooping_integral_species,parameters.values.cols());
	hooping_integral.setZero();
	Eigen::VectorXd v(hooping_integral_species),p(hooping_integral_species);//the final value of hooping integrals between atom1 and atom2.
	v.setZero();
	p.setZero();
	//calculate the group of hooping integrals parameters to be used. 
	int index_of_hooping_parameter_group=2*atoms.element_species+(atom1.element.index*atoms.element_species+atom2.element.index)*hooping_integral_species;
	int index_of_reverse_hooping_parameter_group=2*atoms.element_species+(atom2.element.index*atoms.element_species+atom1.element.index)*hooping_integral_species;
	//calculate the hooping integrals v. 
	for (int i=0;i<hooping_integral_species;i++)
	{
		hooping_integral.row(i)=parameters.values.row(i+index_of_hooping_parameter_group);//read the hooping parameter. 
		double v0=hooping_integral(i,0);
		double r0=hooping_integral(i,1);
		double rc=hooping_integral(i,2);
		double n=hooping_integral(i,3);
		double nc=hooping_integral(i,4);
		v(i)=v0*pow(r0/r,n)*exp(n*(-pow(r/rc,nc)+pow(r0/rc,nc)));;
	}
	
	for (int i=0;i<hooping_integral_species;i++)
	{
		hooping_integral.row(i)=parameters.values.row(i+index_of_reverse_hooping_parameter_group);//read the hooping parameter. 
		double v0=hooping_integral(i,0);
		double r0=hooping_integral(i,1);
		double rc=hooping_integral(i,2);
		double n=hooping_integral(i,3);
		double nc=hooping_integral(i,4);
		p(i)=v0*pow(r0/r,n)*exp(n*(-pow(r/rc,nc)+pow(r0/rc,nc)));
	}

	//calculate the hamiltonian block.
	Eigen::MatrixXcd hamiltonian_block(block_size,block_size);
	hamiltonian_block.setZero();
	hamiltonian_block(0,0)=v(0);//s-s
	hamiltonian_block(0,1)=Azimuth_cos(0)*v(1);//s-px
	hamiltonian_block(0,2)=Azimuth_cos(1)*v(1);//s-py
	hamiltonian_block(0,3)=Azimuth_cos(2)*v(1);//s-pz
	hamiltonian_block(1,0)=-Azimuth_cos(0)*p(1);//px-s
	hamiltonian_block(1,1)=pow(Azimuth_cos(0),2)*v(2)+(1.0-pow(Azimuth_cos(0),2))*v(3);//px-px
	hamiltonian_block(1,2)=Azimuth_cos(0)*Azimuth_cos(1)*(v(2)-v(3));//px-py
	hamiltonian_block(1,3)=Azimuth_cos(0)*Azimuth_cos(2)*(v(2)-v(3));//px-pz
	hamiltonian_block(2,0)=-Azimuth_cos(1)*p(1);//py-s
	hamiltonian_block(2,1)=Azimuth_cos(0)*Azimuth_cos(1)*(p(2)-p(3));//py-px
	hamiltonian_block(2,2)=pow(Azimuth_cos(1),2)*v(2)+(1.0-pow(Azimuth_cos(1),2))*v(3);//py-py
	hamiltonian_block(2,3)=Azimuth_cos(1)*Azimuth_cos(2)*(v(2)-v(3));//py-pz
	hamiltonian_block(3,0)=-Azimuth_cos(2)*p(1);//pz-s
	hamiltonian_block(3,1)=Azimuth_cos(0)*Azimuth_cos(2)*(p(2)-p(3));//pz-px
	hamiltonian_block(3,2)=Azimuth_cos(1)*Azimuth_cos(2)*(p(2)-p(3));//pz-py
	hamiltonian_block(3,3)=pow(Azimuth_cos(2),2)*v(2)+(1.0-pow(Azimuth_cos(2),2))*v(3);//pz-pz

	//calculate the phase factor of bloch wave. 
	double kr=k_points.k_points_coordinates.row(_k_points_index).dot(delta_coordinates);	
	double real=cos(kr),imag=sin(kr);
	std::complex<double> phase_factor(real,imag);

	//hamiltonian block = phase factor * hooping integrals. 
	hamiltonian_block*=phase_factor;
	
	return hamiltonian_block;
}