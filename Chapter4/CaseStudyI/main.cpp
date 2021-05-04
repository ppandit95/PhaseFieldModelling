/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////             CASE STUDY I : Simulation of the Spinodal Decomposition of a binary alloy                   ////////
////////                with explicit Euler Finite Difference Algorithm                                          ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<time.h>
#include<fstream>
#include<cmath>
#include<string>
#include<filesystem>
#include <unistd.h>

#include"src/structures.hpp"
#include"src/SpinodalDecomposition.hpp"


int main(int argc, char **argv) {
	//Get initial wall time
	struct timespec begin,end;
	clock_gettime(CLOCK_REALTIME,&begin);
	SpinodalDecomposition sd;


	//Initial Concentration Profile
	sd.Initial_Profile();


	//Write the initial Profile in a dat file
	sd.write_initial_profile();

	//Evolve with Cahn Hilliard Formulation
	sd.solve();


	clock_gettime(CLOCK_REALTIME,&end);
	std::cout<<"Elapsed Time : "<<(end.tv_sec-begin.tv_sec)<<" seconds "<<(end.tv_nsec - begin.tv_nsec)<<" nanoseconds"<<std::endl;
	return 0;
}
