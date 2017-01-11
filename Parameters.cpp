#include "Parameters.h"
Parameters::Parameters()
{

}

Parameters::~Parameters()
{

}

int Parameters::output( char * filename )
{
	std::ofstream new_param(filename);
	for (int i=0;i<values.rows();i++)
	{
		for (int j=0;j<values.cols();j++)
		{
			new_param<<values(i,j)<<" ";
		}
		for (int j=0;j<fitting_control.cols();j++)
		{
			new_param<<fitting_control(i,j)<<" ";
		}
		if (i!=values.rows()-1)
		{
			new_param<<std::endl;
		}
	}
	new_param.close();
	return 0;
}

