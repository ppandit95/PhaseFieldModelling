/*
 * transient.cpp
 *
 *  Created on: 17-Apr-2021
 *      Author: pushkar
 */
#include"transient.hpp"
#include"structures.hpp"
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



