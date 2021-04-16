///////////////////////////////////////////////////////////////////////////////////
////       1D Transient Heat Conduction with Finite Difference Methods        /////
///////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<fstream>
#include<string>
#include "structures.hpp"
int main(int argc, char **argv) {
//Simulation Parameters
	Params parameters;
	parameters.N = 128;
	parameters.dx = 1.0;
	parameters.mu = 1.0;
	parameters.dt = 0.2;
	parameters.Ti = 1.0;
	parameters.T = new double [parameters.N];
	parameters.temp = new double [parameters.N];

//Setting up initial profile
for(unsigned int i=0;i<parameters.N;i++){
	if(i>((parameters.N/2)-20) && i<((parameters.N/2)+20))
		parameters.T[i] = parameters.Ti;
	else
		parameters.T[i] = 0.0;
}

//Writing the initial Profile
std::ofstream	write_op("Initial_Profile.dat");
for(unsigned int i=0;i<parameters.N;i++)
	write_op << parameters.T[i] <<std::endl;
write_op.close();

for(unsigned int t=0;t<801;t++){
	//Setting up Periodic Boundary
	parameters.temp[0] = parameters.T[0] + parameters.mu*parameters.dt*(1/(parameters.dx*parameters.dx))*(parameters.T[1] - 2.0*parameters.T[0] + parameters.T[parameters.N-1]);
	parameters.temp[parameters.N-1] = parameters.T[parameters.N-1] + parameters.mu*parameters.dt*(1/(parameters.dx*parameters.dx))*(parameters.T[0] - 2.0*parameters.T[parameters.N-1] + parameters.T[parameters.N-2]);

	//Evolving with Jacobi Method
	for(unsigned int i=1;i<parameters.N-1;i++){
		parameters.temp[i] = parameters.T[i] + parameters.mu*parameters.dt*(1/(parameters.dx*parameters.dx))*(parameters.T[i+1] - 2.0*parameters.T[i] + parameters.T[i-1]);
	}

	//Writing Output of parameters.temp values
	if(t%100 == 0){
		std::string filename ="Output-"+std::to_string(t)+".dat";
		std::ofstream out(filename);
		for(unsigned int i=0;i<parameters.N;i++)
			out<<parameters.temp[i]<<std::endl;
		out.close();
	}

	//Retransfering parameters.temp array to T array
	for(unsigned int i=0;i<parameters.N;i++)
		parameters.T[i] = parameters.temp[i];
}

delete[] parameters.T;
delete[] parameters.temp;
	return 0;
}
