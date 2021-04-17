/*
 * transient.cpp
 *
 *  Created on: 17-Apr-2021
 *      Author: pushkar
 */
#include"transient.hpp"
#include"structures.hpp"
#include<fstream>
#include<iostream>
#include<cassert>
TransientProblem::TransientProblem(){
	parameters.N = 128;
	parameters.dx = 1.0;
	parameters.mu = 1.0;
	parameters.dt = 0.2;
	parameters.Ti = 1.0;
	parameters.T = new double [parameters.N];
	parameters.temp = new double [parameters.N];
}
TransientProblem::~TransientProblem(){
	delete[] parameters.T;
	delete[] parameters.temp;
}
Params TransientProblem::GetParameters() const{
	return parameters;
}
void TransientProblem::initialize_field(){
	for(unsigned int i=0;i<parameters.N;i++){
		if(i>((parameters.N/2) - 20) && i<((parameters.N/2) + 20))
			parameters.T[i] = parameters.Ti;
		else
			parameters.T[i] = 0.0;
	}
}
void TransientProblem::Output() const{
	std::ofstream write_output("Initial_Profile.dat");
	assert(write_output.is_open());
	write_output.precision(4);
	for(unsigned int i=0;i<parameters.N;i++)
		write_output<<parameters.T[i]<<std::endl;
	write_output.close();
}
void TransientProblem::setup_periodic_BC(){
	parameters.temp[0] = parameters.T[0] + parameters.mu*parameters.dt*(1/(parameters.dx*parameters.dx))*
						(parameters.T[1] - 2.0*parameters.T[0] + parameters.T[parameters.N-1]);

	parameters.temp[parameters.N-1] = parameters.T[parameters.N-1] + parameters.mu*parameters.dt*
			                          (1/(parameters.dx*parameters.dx))*(parameters.T[0] - 2.0*parameters.T[parameters.N-1]
									+ parameters.T[parameters.N-2]);

}
void TransientProblem::Jacobi_iteration(){
	for(unsigned int i=1;i<parameters.N-1;i++){
			parameters.temp[i] = parameters.T[i] + parameters.mu*parameters.dt*(1/(parameters.dx*parameters.dx))
					*(parameters.T[i+1] - 2.0*parameters.T[i] + parameters.T[i-1]);
		}
}
void TransientProblem::Output(unsigned int t) const{
	if(t%100 == 0){
			std::string filename ="Output-"+std::to_string(t)+".dat";
			std::ofstream out(filename);
			for(unsigned int i=0;i<parameters.N;i++)
				out<<parameters.temp[i]<<std::endl;
			out.close();
		}
}
void TransientProblem::retransfer(){
	for(unsigned int i=0;i<parameters.N;i++)
			parameters.T[i] = parameters.temp[i];
}
void TransientProblem::solve(){
	for(unsigned int t=0;t<801;t++){
		//Setting up Periodic Boundary
		setup_periodic_BC();
		//Evolving with Jacobi Method
		Jacobi_iteration();

		//Writing Output of prob.GetParameters().temp values
		Output(t);

		//Retransfering prob.GetParameters().temp array to T array
		retransfer();
	}
}


