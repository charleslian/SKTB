
#include <iostream>
#include "Fitting.h"
Fitting::Fitting()
{

}

Fitting::~Fitting()
{

}

int Fitting::fitting(double criteria) 
{
	create_simplex();

	int N=simplex.cols();
	Eigen::VectorXd function_values(N+1);function_values.setZero();
	Eigen::VectorXd function_values_temp(N+1);function_values_temp.setZero();
	Eigen::RowVectorXd center(N);center.setZero();
	Eigen::RowVectorXd try_simplex_vertex(N);try_simplex_vertex.setZero();

	for (int k=0;k<N+1;k++)
	{
		center+=simplex.row(k)/(N+1);
		simplex_to_parameters(simplex.row(k));
		function_values(k)=fitting_function();
	}
	while (true)
	{
		//find the maximum, minimum, next-maximum values and corresponding vertex.
		int minimum_index=0;
		int maximum_index=0;
		int next_maximum_index=0;
		
		double minimum_value=function_values.minCoeff(&minimum_index);
		double maximum_value=function_values.maxCoeff(&maximum_index);
		
		//set the maximum value of function_values_temp to minimum_value. 
		//find maximum of function_values_temp to get the next-maximum value.
		function_values_temp=function_values;function_values_temp(maximum_index)=minimum_value;
		double next_maximum_value=function_values_temp.maxCoeff(&next_maximum_index);

		
		//calculate the difference between maximum and minimum. 
		double difference=2.0*fabs(maximum_value-minimum_value)/(fabs(maximum_value)+fabs(maximum_value));
		std::cout<<" minimum value:"<<minimum_value<<" maximum value:"<<maximum_value<<" next maximum value:"<<next_maximum_value<<" tolerance:"<<difference<<std::endl;
		
		//if difference is smaller than the criteria, then stop. 
		if (difference<criteria)
		{
			simplex_to_parameters(simplex.row(minimum_index));
			break;
		}

		//try to reflect the highest point. 
		try_simplex_vertex=2.0*center-(2.0/(N+1)+1)*simplex.row(maximum_index);
		simplex_to_parameters(try_simplex_vertex);
		double try_value=fitting_function();

		//if (try_value!=try_value)
		//{
		//	try_value=1000;
		//}
		//std::cout<<try_value<<std::endl;
		//if reflection value is better than the highest point then replace the highest point.  
		if (try_value<function_values(maximum_index))
		{
			center+=(try_simplex_vertex-simplex.row(maximum_index))/(N+1);
			simplex.row(maximum_index)=try_simplex_vertex;
			function_values(maximum_index)=try_value;
		}

		//if reflection value is even better than the lowest point.  
		if (try_value<function_values(minimum_index))
		{
			//try to expand along the direction.  
			try_simplex_vertex=-1.0*center+(1.0/(N+1)+2.0)*simplex.row(maximum_index);
			simplex_to_parameters(try_simplex_vertex);
			try_value=fitting_function();

			//if the expansion value is better than the reflection point
			//replace the highest point with the expansion value.
			if (try_value<function_values(maximum_index))
			{
				center+=(try_simplex_vertex-simplex.row(maximum_index))/(N+1);
				simplex.row(maximum_index)=try_simplex_vertex;
				function_values(maximum_index)=try_value;
			}
		}

		//if reflection value is worse than the next highest point.  
		else if(try_value>function_values(next_maximum_index))
		{
			//try to contract along the direction. 
			try_simplex_vertex=0.5*center+(0.5/(N+1)-0.5)*simplex.row(maximum_index);
			simplex_to_parameters(try_simplex_vertex);
			try_value=fitting_function();
			//if the contraction value is worse than the highest point. 
			if (try_value>function_values(maximum_index))
			{
				//contract around the lowest point. 
				for (int i=0;i<N+1;i++)
				{
					if (i!=minimum_index)
					{
						try_simplex_vertex=0.5*(simplex.row(i)+simplex.row(minimum_index));
						center+=(try_simplex_vertex-simplex.row(i))/(N+1);
						simplex.row(i)=try_simplex_vertex;
						simplex_to_parameters(try_simplex_vertex);
						function_values(i)=fitting_function();
					}
				}
			}
			//if the contraction value is better than the highest point, then replace the highest point.
			else
			{
				center+=(try_simplex_vertex-simplex.row(maximum_index))/(N+1);
				simplex.row(maximum_index)=try_simplex_vertex;
				function_values(maximum_index)=try_value;
			}
		}
// 		if (try_value<function_values(maximum_index))
// 		{
// 			center+=(try_simplex_vertex-simplex.row(maximum_index))/(N+1);
// 			simplex.row(maximum_index)=try_simplex_vertex;
// 		}	
	}
	return 0;
}

int Fitting::create_simplex() 
{
	//set the simplex from simplex without control.
	//if the fitting_control(i,j)==0, then omit the parameter (i,j).
	int N=0;
	//count the number of 1 in fitting_control.
	for (int i=0;i<parameters.values.rows();i++)
	{
		for (int j=0;j<parameters.values.cols();j++)
		{
			if (parameters.fitting_control(i,j)!=0)
			{
				N++;
			}
		}
	}

	//read the parameters at the position (i,j) where fitting_control(i,j)==0.
	simplex=Eigen::MatrixXd(N+1,N);
	for (int k=0;k<N+1;k++)
	{
		int m=0;
		Eigen::MatrixXd simplex_vertex_without_control=set_simplex_vertex_without_control(k);
		for (int i=0;i<parameters.values.rows();i++)
		{
			for (int j=0;j<parameters.values.cols();j++)
			{				
				if (parameters.fitting_control(i,j)!=0)
				{					
					simplex(k,m)=simplex_vertex_without_control(i,j);
					m++;
				}
			}
		}
	}
	//std::cout<<simplex<<std::endl;
	return 0;
}

int Fitting::simplex_to_parameters(Eigen::RowVectorXd try_simplex_vertex) 
{
	//generate parameters from simplex.
	//if the fitting_control(i,j)==0, then do not change the parameter (i,j) and go to next.
	int k=0;
	for (int i=0;i<parameters.values.rows();i++)
	{
		for (int j=0;j<parameters.values.cols();j++)
		{
			if (parameters.fitting_control(i,j)!=0)
			{
				parameters.values(i,j)=try_simplex_vertex(k);
				k++;
			}
		}
	}
	return 0;
}

Eigen::MatrixXd Fitting::set_simplex_vertex_without_control(int _index_of_vertex) 
{
	int NX=parameters.values.rows();
	int NY=parameters.values.cols();
	int N=NX*NY;
	Eigen::MatrixXd simplex_vertex_without_control(NX,NY);
	return simplex_vertex_without_control;
}

double Fitting::fitting_function()
{
	return 0;
}
