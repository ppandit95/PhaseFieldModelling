///////////////////////////////////////////////////////////////////////////////////
////       1D Transient Heat Conduction with Finite Difference Methods        /////
///////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<fstream>
#include<string>
int main(int argc, char **argv) {
//Simulation Parameters

unsigned int N = 128;
double dx = 1;
double mu = 1.0;
double dt = 0.2;
double Ti = 1.0;
double T[N],temp[N];

//Setting up initial profile
for(unsigned int i=0;i<N;i++){
	if(i>((N/2)-20) && i<((N/2)+20))
		T[i] = Ti;
	else
		T[i] = 0.0;
}

//Writing the initial Profile
std::ofstream	write_op("Initial_Profile.dat");
for(unsigned int i=0;i<N;i++)
	write_op << T[i] <<std::endl;
write_op.close();

for(unsigned int t=0;t<801;t++){
	//Setting up Periodic Boundary
	temp[0] = T[0] + mu*dt*(1/(dx*dx))*(T[1] - 2.0*T[0] + T[N-1]);
	temp[N-1] = T[N-1] + mu*dt*(1/(dx*dx))*(T[0] - 2.0*T[N-1] + T[N-2]);

	//Evolving with Jacobi Method
	for(unsigned int i=1;i<N-1;i++){
		temp[i] = T[i] + mu*dt*(1/(dx*dx))*(T[i+1] - 2.0*T[i] + T[i-1]);
	}

	//Writing Output of temp values
	if(t%100 == 0){
		std::string filename ="Output-"+std::to_string(t)+".dat";
		std::ofstream out(filename);
		for(unsigned int i=0;i<N;i++)
			out<<temp[i]<<std::endl;
		out.close();
	}

	//Retransfering temp array to T array
	for(unsigned int i=0;i<N;i++)
		T[i] = temp[i];
}

	return 0;
}
