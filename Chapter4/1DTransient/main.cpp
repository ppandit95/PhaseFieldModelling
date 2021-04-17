///////////////////////////////////////////////////////////////////////////////////
////       1D Transient Heat Conduction with Finite Difference Methods        /////
///////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<fstream>
#include<string>
#include "structures.hpp"
#include"transient.hpp"
int main(int argc, char **argv) {
	TransientProblem prob;


//Setting up initial profile
for(unsigned int i=0;i<prob.GetParameters().N;i++){
	if(i>((prob.GetParameters().N/2)-20) && i<((prob.GetParameters().N/2)+20))
		prob.GetParameters().T[i] = prob.GetParameters().Ti;
	else
		prob.GetParameters().T[i] = 0.0;
}

//Writing the initial Profile
std::ofstream	write_op("Initial_Profile.dat");
for(unsigned int i=0;i<prob.GetParameters().N;i++)
	write_op << prob.GetParameters().T[i] <<std::endl;
write_op.close();

for(unsigned int t=0;t<801;t++){
	//Setting up Periodic Boundary
	prob.GetParameters().temp[0] = prob.GetParameters().T[0] + prob.GetParameters().mu*prob.GetParameters().dt*(1/(prob.GetParameters().dx*prob.GetParameters().dx))*(prob.GetParameters().T[1] - 2.0*prob.GetParameters().T[0] + prob.GetParameters().T[prob.GetParameters().N-1]);
	prob.GetParameters().temp[prob.GetParameters().N-1] = prob.GetParameters().T[prob.GetParameters().N-1] + prob.GetParameters().mu*prob.GetParameters().dt*(1/(prob.GetParameters().dx*prob.GetParameters().dx))*(prob.GetParameters().T[0] - 2.0*prob.GetParameters().T[prob.GetParameters().N-1] + prob.GetParameters().T[prob.GetParameters().N-2]);

	//Evolving with Jacobi Method
	for(unsigned int i=1;i<prob.GetParameters().N-1;i++){
		prob.GetParameters().temp[i] = prob.GetParameters().T[i] + prob.GetParameters().mu*prob.GetParameters().dt*(1/(prob.GetParameters().dx*prob.GetParameters().dx))*(prob.GetParameters().T[i+1] - 2.0*prob.GetParameters().T[i] + prob.GetParameters().T[i-1]);
	}

	//Writing Output of prob.GetParameters().temp values
	if(t%100 == 0){
		std::string filename ="Output-"+std::to_string(t)+".dat";
		std::ofstream out(filename);
		for(unsigned int i=0;i<prob.GetParameters().N;i++)
			out<<prob.GetParameters().temp[i]<<std::endl;
		out.close();
	}

	//Retransfering prob.GetParameters().temp array to T array
	for(unsigned int i=0;i<prob.GetParameters().N;i++)
		prob.GetParameters().T[i] = prob.GetParameters().temp[i];
}


	return 0;
}
