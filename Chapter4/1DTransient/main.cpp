///////////////////////////////////////////////////////////////////////////////////
////       1D Transient Heat Conduction with Finite Difference Methods        /////
///////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<fstream>
#include<string>
#include "src/structures.hpp"
#include"src/transient.hpp"
int main(int argc, char **argv) {
	TransientProblem prob;


	//Setting up initial profile
	prob.initialize_field();

	//Writing the initial Profile
	prob.Output();

	prob.solve();


	return 0;
}
