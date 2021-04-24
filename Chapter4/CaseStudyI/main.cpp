/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////             CASE STUDY I : Simulation of the Spinodal Decomposition of a binary alloy                   ////////
////////                with explicit Euler Finite Difference Algorithm                                          ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<time.h>
#include<fstream>
int main(int argc, char **argv) {
	//Get initial wall time
	struct timespec begin,end;
	clock_gettime(CLOCK_REALTIME,&begin);

	//Open output file to write energy functional values at each timestep
	//std::ofstream out("time_energ.dat");

	//Simulation cell parameters
	int Nx = 64;
	int Ny = 64;
	double dx = 1.0;
	double dy = 1.0;
	double** con,mu;

	//Time Integration Parameters
	int nstep = 10000;
	int nprint = 50;
	double dtime = 0.01;
	double ttime = 0.0;

	//Material Specific Parameters
	double c0 = 0.4;
	double mobility = 1.0;
	double grad_coef = 0.5;

	//Initial concentration profile
	srand( (unsigned)time( NULL ) );
	con = new double*[Nx];
	mu = new double* [Nx];
	for(int i=0;i<Nx;i++){
		con[i] = new double[Ny];
		mu[i] = new double[Ny];
	}
	double noise = 0.02;
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			con[i][j] = c0 + noise*(0.5 - (float) rand()/RAND_MAX);
		}
	}

	//Write the initial Profile in a dat file
	std::ofstream initProf("Initial_Profile.dat");
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++)
			initProf << con[i][j] << "\t";
		initProf << std::endl;
	}

	//Evolve with Cahn Hilliard Formulation
	for(int istep=0;istep<nstep;istep++){
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				//Compute the chemical potential at each time step and each point

			}
		}



	}
	for(int i=0;i<Nx;i++){
		delete[] con[i];
		delete[] mu[i];
	}
	delete[] mu;
	delete[] con;
	return 0;
}
