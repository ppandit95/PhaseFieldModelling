///////////////////////////////////////////////////////////////////////////////////
////       1D Transient Heat Conduction with Finite Difference Methods        /////
///////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<fstream>
int main(int argc, char **argv) {
//Simulation Parameters

unsigned int N = 128;
double dx = 1;
double mu = 1.0;
double dt = 0.2;
double Ti = 1.0;
double T[N];

//Setting up initial profile
for(unsigned int i=0;i<N;i++){
	if(i>((N/2)-40) && i<((N/2)+40))
		T[i] = Ti;
	else
		T[i] = 0.0;
}

//Writing the initial Profile
std::ofstream	write_op("Initial_Profile.dat");
for(unsigned int i=0;i<N;i++)
	write_op << T[i] <<std::endl;
write_op.close();

	return 0;
}
